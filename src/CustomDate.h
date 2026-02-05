#pragma once
#include <ctime>
#include <string>

class CustomDate {
   public:
    static const double matlabd_epoch;
    CustomDate() {}
    CustomDate(double matlabd);
    CustomDate(const CustomDate& rhs);
    CustomDate& operator=(const CustomDate& rhs);
    CustomDate operator+(uint32_t delta_seconds);
    CustomDate operator-(uint32_t delta_seconds);
    CustomDate add_month() const;

    std::string to_date_string() const;
    std::string to_string() const;
    CustomDate to_bom() const;  // to beginning of month
    CustomDate to_eom() const;  // to end of month
    bool operator==(const CustomDate& rhs);
    bool operator!=(const CustomDate& rhs);
    bool operator<(const CustomDate& rhs);
    bool operator<=(const CustomDate& rhs);
    bool operator>(const CustomDate& rhs);
    bool operator>=(const CustomDate& rhs);
   private:
    std::tm datetime;
    std::time_t unix_seconds = 0;
    static const int max_days[12];
    bool is_leap_year() const;
    int eom() const;
};
