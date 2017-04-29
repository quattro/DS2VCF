#include <algorithm>
#include <cmath>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>


#include "gzstream.h"
#include "readerwriterqueue.h"
#include "atomicops.h"

namespace mc = moodycamel;

template<typename Out> void split(const std::string &s, char delim, Out result);
std::vector<std::string> split(const std::string &s, char delim);
bool startsWith(const std::string& haystack, const std::string& needle);
double PLtoPP(std::string x);


int main (int argc, char *argv[]) {

    const int REF = 3;
    const int ALT = 4;
    const int FORMAT = 8;
    const int ENTRY_START = FORMAT + 1;
    const size_t Q_SIZE = 1000; 
    const std::string QUIT_STR = "~";


    igzstream infile(argv[1]);
    ogzstream outfile(argv[2]);

    mc::BlockingReaderWriterQueue<std::string> lines(Q_SIZE);


    // set up consumer thread
    std::thread reader([&]() {
        double pls [3];
        std::string line;
        std::string entry;
        std::stringstream ss_line;

        ss_line.setf(std::ios::fixed,std::ios::floatfield);
        ss_line.precision(2);

        while (true) {
            lines.wait_dequeue(line);

            if (line.compare(QUIT_STR) == 0) {
                // output any remaining strings
                outfile << ss_line.rdbuf();
                break;
            }

            if (startsWith(line, "#")) {
                ss_line << line << std::endl;
                continue;
            }

            std::vector<std::string> row = split(line, '\t');

            if (row[REF].length() > 1 || row[ALT].length() > 1) {
                continue;
            }
            for (int i = 0; i < FORMAT; i++) {
                ss_line << row[i] << '\t';
            }
            ss_line << row[FORMAT] + ":DS" << '\t';

            for (int i = ENTRY_START; i < row.size() - 1; i++) {
                std::string token = row[i].substr(row[i].find_last_of(':') + 1);
                std::stringstream ss;
                ss.str(token);
    
                for (int j = 0; j < 3; j++) {
                    std::getline(ss, entry, ',');
                    pls[j] = PLtoPP(entry);
                }
    
                double dose = (pls[1] + 2 * pls[2]) / (pls[0] + pls[1] + pls[2]);
                ss_line << row[i] + ":" << dose << '\t';
            }
            int i = row.size() - 1;
            std::string token = row[i].substr(row[i].find_last_of(':') + 1);
            std::stringstream ss;
            ss.str(token);
    
            for (int j = 0; j < 3; j++) {
                std::getline(ss, entry, ',');
                pls[j] = PLtoPP(entry);
            }
    
            double dose = (pls[1] + 2 * pls[2]) / (pls[0] + pls[1] + pls[2]);
            ss_line << row[i] + ":" << dose << std::endl;

            outfile << ss_line.rdbuf();
            ss_line.str(std::string());
        }
    });

    std::string line;
    while (std::getline(infile, line)) {
        while(!lines.try_enqueue(line)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    }
    while(!lines.try_enqueue(QUIT_STR)) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }

    reader.join();
    outfile.close();
    infile.close();

    return 0;
}

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

bool startsWith(const std::string& haystack, const std::string& needle) {
        return needle.length() <= haystack.length() 
                    && std::equal(needle.begin(), needle.end(), haystack.begin());
}


double PLtoPP(std::string x) {
    std::string::size_type sz; 
    return std::pow(10.0, std::stod(x, &sz) / -10.0);
}
