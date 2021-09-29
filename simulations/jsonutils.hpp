#include <iostream>
#include <iomanip>

//Function to write jsonArray from any iterable to any stream (file or cout)
template <typename T, typename Stream> void  jsonArray(T& toPrint, Stream &stream, long double Norm=1)
{
    stream<<"[";
    bool firstVal=true;
    for( auto val: toPrint)
    {
        if (!firstVal)
            stream<<",";
        firstVal=false;

        if(Norm == 1)
            stream << std::setprecision(15) << val;
        else
            stream << std::setprecision(15) << val / Norm;

    }
    stream<<"]";
    stream.flush();
}  

//Function to write jsonArrayofArrays from any iterable of iterables to any stream (file or cout)
template <typename T, typename Stream> void  jsonArrayofArrays(T& toPrint, Stream &stream, long double Norm=1)
{
    stream<<"[";
    bool firstVal=true;
    for( auto val: toPrint)
    {
        if (!firstVal)
            stream<<",";
        firstVal=false;
        jsonArray(val,stream,Norm);
        
    }
    stream<<"]";
    stream.flush();
}  


//Function to write jsonArray of [first,second] from any iterable of std::pairs to any stream (file or cout)
template <typename T, typename Stream> void  jsonPairArray(T toPrint, Stream &stream, long double Norm=1)
{
    stream<<"[";
    bool firstVal=true;
    for( auto val: toPrint)
    {
        if (!firstVal)
            stream<<",";
        firstVal=false;

        if(Norm == 1)
            stream << "[" << std::setprecision(15) << val.first << "," << val.second <<"]";
        else
            stream << "[" << std::setprecision(15) << val.first/Norm << "," << val.second/Norm <<"]";

    }
    stream<<"]";
    stream.flush();
}  


//Function to write json dictionary from any map to any stream (file or cout)
template <typename T, typename Stream> void  jsonMap(T toPrint, Stream &stream, long double Norm=1)
{
    stream<<"{";
    bool firstVal=true;
    for( auto val: toPrint)
    {
        if (!firstVal)
            stream<<",";
        firstVal=false;
        stream << "\"" << val.first << "\" : ";
        if(Norm == 1)
            stream << std::setprecision(15) << val.second;
        else
            stream <<std::setprecision(15) << val.second / Norm;

    }
    stream<<"}";
    stream.flush();
}  

//Function to write json dictionary of arrays from any map of arrays to any stream (file or cout)
template <typename T, typename Stream> void  jsonMapOfArrays(T toPrint, Stream &stream, long double Norm=1)
{
    stream<<"{";
    bool firstVal=true;
    for( auto val: toPrint)
    {
        if (!firstVal)
            stream<<",";
        firstVal=false;
        stream << "\"" << val.first << "\" : ";
        jsonArray(val.second,stream);
    }
    stream<<"}";
    stream.flush();
}  