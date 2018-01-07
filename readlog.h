#ifndef READLOG_H
#define READLOG_H
#include <string>
#include "traiettoria.h"
#include <vector>


class ReadLog
{
public:
    ReadLog(std::string filename,Traiettoria * t=0,unsigned int skip=1);
    ~ReadLog();
    double * line(unsigned int index);
    unsigned int n_timestep(){return data.size()/data_size;}
    unsigned int n_data(){return headers.size()-1;}
    unsigned int timestep(unsigned int index);
    unsigned int get_natoms(){return 1728;}
    std::pair<unsigned int,bool> get_index_of(std::string header);
private:
    Traiettoria * traiettoria;
    std::vector<std::string> headers;
    std::vector<double > data;
    std::vector<unsigned int > timesteps;
    unsigned int skip,size,step_index;
    size_t data_size;
    bool if_only_numbers(std::string str);
    unsigned int natoms;

};

#endif // READLOG_H
