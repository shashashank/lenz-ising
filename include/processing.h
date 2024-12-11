// #include<cstdlib>
// #include<stdexcept>
// #include<cstdio>
// #include<cmath>
// #include<iostream>
// #include<fstream>
// #include<vector>
// #include<array>
// #include<random>
// #include<dirent.h>
// #include<string> //for string
// #include<sstream> //covert into to string
// #include<sys/types.h> //mkdir
// #include<sys/stat.h>   //mkdir
// #include<iomanip> // setprecision
// #include<omp.h>
// #include<list>
// #include<algorithm>
// #include <cstddef>
// #include <iterator>
// #include<chrono>
// #include <iomanip>
// #include <ctime>
// #include <filesystem>


// void extractEandM(std::ifstream& infile, std::vector<double>& e, std::vector<double>& m);
// double autoCorrelate(const std::vector<double>& obsVec, std::vector<double>& acVec, double kMax);
// double binningError(const std::vector<double>& obsVec, const size_t nB);
// double jackknifeError(const std::vector<double>& obsVec, const size_t nB);
// void binAnalysisFromFile(std::fstream& infile, const size_t obsColumn);
// double binningErrorFromFile(std::fstream&infile, std::vector<double>& binVec, const size_t obsColumn, const size_t nB);
// void errorToFile(std::ifstream& indata, std::ofstream& odata, std::vector<double>& eVec, std::vector<double>& mVec);

// inline double meanOfVector(std::vector<double>::const_iterator start, const std::vector<double>::const_iterator end){
//     const size_t n = end - start;
//     double mean = std::accumulate(start, end, 0.0);
//     return mean / n;
// }

// inline double meanOfVector(const std::vector<double> & obsVec){
//     double mean = 0.0;

//     for (size_t i = 0; i < obsVec.size(); ++i)
//     {
//         mean += obsVec[i];
//     }

//     return mean / obsVec.size();
// }

// inline double varianceOfVector(const std::vector<double>& obsVec){
//     const double mean = meanOfVector(obsVec);
//     double variance = 0.0, diff = 0.0;

//     for (size_t i = 0; i < obsVec.size(); ++i)
//     {
//         diff = obsVec[i] - mean;
//         variance += diff * diff;
//     }

//     return variance/(obsVec.size());
// }

// // inline double varianceOfVector(const std::vector<double>& obsVec){
// //     double c = 0, b = 0, d=0;
// //     size_t N = obsVec.size();
// //     auto oPtr = obsVec.begin();
// //     while(oPtr != obsVec.end()){
// //         d = *oPtr++;
// //         c += d;
// //         b += d * d;
// //     }
// //     c = c*c/(N*N);
// //     b = b/N;
// //     return b - c;
// // }

// inline size_t mod(size_t a, size_t b)
// {
//     int ret = a % b;
//     if (ret<0)
//         ret+=b;
//     return ret;
// }

// inline int readFile(const std::string fileName, std::vector<std::string>& data){
//     std::ifstream infile(fileName);
//     std::string line;
//     while (getline(infile, line)){
//         data.push_back(line);
//     }
//     return data.size();
// }

// inline int readFile(std::ifstream& file, std::vector<std::string>& data){
//     std::string line;
//     while (getline(file, line)){
//         data.push_back(line);
//     }
//     return data.size();
// }

// class progress_bar
// {
//     static const auto overhead = sizeof " [100%]";
//     std::ostream& os;
//     const std::size_t bar_width;
//     std::string message;
//     const std::string full_bar;

//  public:
//     progress_bar(std::ostream& os, std::size_t line_width,
//                  std::string message_, const char symbol = '.')
//         : os{os},
//           bar_width{line_width - overhead},
//           message{std::move(message_)},
//           full_bar{std::string(bar_width, symbol) + std::string(bar_width, ' ')}
//     {
//         if (message.size()+1 >= bar_width || message.find('\n') != message.npos) {
//             os << message << '\n';
//             message.clear();
//         } else {
//             message += ' ';
//         }
//         write(0.0);
//     }

//     // not copyable
//     progress_bar(const progress_bar&) = delete;
//     progress_bar& operator=(const progress_bar&) = delete;

//     ~progress_bar()
//     {
//         write(1.0);
//         os << '\n';
//     }

//     void write(double fraction);
// };

// inline void progress_bar::write(double fraction)
// {
//     // clamp fraction to valid range [0,1]
//     if (fraction < 0)
//         fraction = 0;
//     else if (fraction > 1)
//         fraction = 1;

//     auto width = bar_width - message.size();
//     auto offset = bar_width - static_cast<unsigned>(width * fraction);

//     os << '\r' << message;
//     os.write(full_bar.data() + offset, width);
//     os << " [" << std::setw(3) << static_cast<int>(100*fraction) << "%] " << std::flush;
// }
