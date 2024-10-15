#include "CustomDate.h"

#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <sstream>
#include <tuple>
#include <iomanip>

#ifdef _MSC_VER
    #define _CRT_SECURE_NO_WARNINGS
    // Define gmtime_r as gmtime_s on MSVC
    #define gmtime_r( ARG1, ARG2 ) gmtime_s( ARG2, ARG1 )
#else
    // On POSIX systems (like Linux), gmtime_r already exists
    #include <time.h>
#endif


const double CustomDate::matlabd_epoch = 719529.0;
const int CustomDate::max_days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

CustomDate::CustomDate(double matlabd) {
    if (matlabd < matlabd_epoch) {
        std::cout << "Input date is smaller than the epoch." << std::endl;
        exit(EXIT_FAILURE);
    }

    unix_seconds = static_cast<time_t>(
        (matlabd - matlabd_epoch) * 24. * 3600. + 0.5);
    gmtime_r(&unix_seconds, &datetime);
}

CustomDate::CustomDate(const CustomDate& rhs) : datetime(rhs.datetime), unix_seconds(rhs.unix_seconds) {}

CustomDate& CustomDate::operator=(const CustomDate& rhs) {
    unix_seconds = rhs.unix_seconds;
    datetime = rhs.datetime;
    return *this;
}

CustomDate CustomDate::operator+(uint32_t delta_seconds)
{
    CustomDate answer;
    answer.unix_seconds = unix_seconds + delta_seconds;
    gmtime_r(&answer.unix_seconds, &answer.datetime);
    return answer;
}
CustomDate CustomDate::operator-(uint32_t delta_seconds)
{
    CustomDate answer;
    answer.unix_seconds = unix_seconds - delta_seconds;
    gmtime_r(&answer.unix_seconds, &answer.datetime);
    return answer;
}

CustomDate CustomDate::add_month() const {
    CustomDate answer;
    int delta_days         = max_days[datetime.tm_mon];
    int days_in_next_month = max_days[(datetime.tm_mon + 1) % 12];
    
    if(this->is_leap_year()){
        if(datetime.tm_mon == 0){
            days_in_next_month += 1;
        }
        else if(datetime.tm_mon == 1){
            delta_days += 1;
        }
    }

    delta_days -= std::max(0, datetime.tm_mday - days_in_next_month);
    answer.unix_seconds = unix_seconds + delta_days * 24 * 3600;
    gmtime_r(&answer.unix_seconds, &answer.datetime);

    return answer;
}

std::string CustomDate::to_date_string() const {
    std::stringstream str_stream;
    str_stream << std::setfill('0')
        << datetime.tm_year + 1900 
        << std::setw(2) << datetime.tm_mon + 1
        << std::setw(2) << datetime.tm_mday;

    return str_stream.str();
}
std::string CustomDate::to_string() const {
    std::stringstream str_stream;
    str_stream << std::setfill('0')
        << datetime.tm_year + 1900 << '-'
        << std::setw(2) << datetime.tm_mon + 1 << '-'
        << std::setw(2) << datetime.tm_mday << 'T'
        << std::setw(2) << datetime.tm_hour << ':'
        << std::setw(2) << datetime.tm_min << ':'
        << std::setw(2) << datetime.tm_sec;

    return str_stream.str();
}

CustomDate CustomDate::to_bom() const {
    CustomDate answer;
    answer.unix_seconds = unix_seconds
        - (datetime.tm_mday - 1) * 24 * 3600
        - datetime.tm_hour * 3600
        - datetime.tm_min * 60
        - datetime.tm_sec;
    gmtime_r(&answer.unix_seconds, &answer.datetime);
    return answer;
}

CustomDate CustomDate::to_eom() const {
    CustomDate answer;
    int new_mday = this->eom();
    answer.unix_seconds = unix_seconds 
        + (new_mday - datetime.tm_mday) * 24 * 3600 // seconds in delta days
        + (23 - datetime.tm_hour) * 3600 // seconds in delta hours
        + (59 - datetime.tm_min) * 60 // seconds in delta minutes
        + (59 - datetime.tm_sec);

    gmtime_r(&answer.unix_seconds, &answer.datetime);
    return answer;
}

bool CustomDate::is_leap_year() const {
    int y = datetime.tm_year + 1900;
    if (((y % 4 == 0) && (y % 100 != 0)) || (y % 400 == 0))
        return true;

    return false;
}

int CustomDate::eom() const {
    int new_day = max_days[datetime.tm_mon];
    if (is_leap_year() && datetime.tm_mon == 1) {
        new_day++;
    }
    return new_day;
}

bool CustomDate::operator==(const CustomDate& rhs) {
    return unix_seconds == rhs.unix_seconds;
}
bool CustomDate::operator!=(const CustomDate& rhs) {
    return unix_seconds != rhs.unix_seconds;
}
bool CustomDate::operator<(const CustomDate& rhs) {
    return unix_seconds < rhs.unix_seconds;
}
bool CustomDate::operator<=(const CustomDate& rhs) {
    return unix_seconds <= rhs.unix_seconds;
}
bool CustomDate::operator>(const CustomDate& rhs) {
    return unix_seconds > rhs.unix_seconds;
}
bool CustomDate::operator>=(const CustomDate& rhs) {
    return unix_seconds >= rhs.unix_seconds;
}
