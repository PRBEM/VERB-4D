#pragma once

#include <string>

class CustomDate {
public:

    static const double matlabd_epoch;

    CustomDate(int y=0, int m=1, int d=1) 
        :year(y), month(m), day(d) {}
    CustomDate(double matlabd);
    CustomDate(const CustomDate& rhs);
    CustomDate& operator=(const CustomDate& rhs);

    int get_second() const {return second;}
    int get_minute() const {return minute;}
    int get_hour() const {return hour;}
    int get_day() const {return day;}
    int get_month() const {return month;}
    int get_year() const {return year;}

    CustomDate add_month(int add = 1) const;

    std::string to_string() const;
    CustomDate to_bom() const; // to beginning of month
    CustomDate to_eom() const; // to end of month
private:
    int year;
    int month;
    int day;
    int hour, minute, second;

    static const int max_days[12];

    static bool is_leap_year(int year);

    int eom() const;
};

bool operator==(const CustomDate& rhs, const CustomDate& lhs);
bool operator!=(const CustomDate& rhs, const CustomDate& lhs);
bool operator<(const CustomDate& rhs, const CustomDate& lhs);
bool operator<=(const CustomDate& rhs, const CustomDate& lhs);
bool operator>(const CustomDate& rhs, const CustomDate& lhs);
bool operator>=(const CustomDate& rhs, const CustomDate& lhs);
