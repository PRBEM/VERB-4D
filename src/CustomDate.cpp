#include "CustomDate.h"

#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <tuple>

const double CustomDate::matlabd_epoch = 719529.0;
const int CustomDate::max_days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

CustomDate::CustomDate(double matlabd) {
    if (matlabd < matlabd_epoch) {
        std::cout << "Input date is smaller than the epoch." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::time_t unix_time = static_cast<time_t>(
        (matlabd - matlabd_epoch) * 24. * 3600. + 0.5);
    std::tm* tm_tmp = gmtime(&unix_time);

    year = 1900 + tm_tmp->tm_year;
    month = 1 + tm_tmp->tm_mon;
    day = tm_tmp->tm_mday;
    hour = tm_tmp->tm_hour;
    minute = tm_tmp->tm_min;
    second = tm_tmp->tm_sec;
}

CustomDate::CustomDate(const CustomDate& rhs) {
    second = rhs.second;
    minute = rhs.minute;
    hour = rhs.hour;
    day = rhs.day;
    month = rhs.month;
    year = rhs.year;
}

CustomDate& CustomDate::operator=(const CustomDate& rhs) {
    second = rhs.second;
    minute = rhs.minute;
    hour = rhs.hour;
    day = rhs.day;
    month = rhs.month;
    year = rhs.year;

    return *this;
}

CustomDate CustomDate::add_month(int m) const {
    CustomDate date_tmp = *this;

    date_tmp.month += m;
    if (date_tmp.month > 12) {
        date_tmp.year += (date_tmp.month / 12);
        date_tmp.month = date_tmp.month % 12;
    } else if (date_tmp.month <= 0) {
        date_tmp.year += (date_tmp.month - 12) / 12;
        date_tmp.month = 12 - abs(date_tmp.month % 12);
    }

    int max_day = max_days[date_tmp.month - 1];
    if (is_leap_year(date_tmp.year) && date_tmp.month == 2) {
        max_day += 1;
    }

    if (date_tmp.day > max_day) {
        date_tmp.day = max_day;
    }

    return date_tmp;
}

std::string CustomDate::to_full_string() const {
    std::string str_rep = std::to_string(year);

    if (month < 10) {
        str_rep += ("0" + std::to_string(month));
    } else {
        str_rep += std::to_string(month);
    }

    if (day < 10) {
        str_rep += ("0" + std::to_string(day));
    } else {
        str_rep += std::to_string(day);
    }

    if (hour < 10) {
        str_rep += (" 0" + std::to_string(hour));
    } else {
        str_rep += " " + std::to_string(hour);
    }

    if (minute < 10) {
        str_rep += (":0" + std::to_string(minute));
    } else {
        str_rep += ":" + std::to_string(minute);
    }

    if (second < 10) {
        str_rep += (":0" + std::to_string(second));
    } else {
        str_rep += ":" + std::to_string(second);
    }

    return str_rep;
}

std::string CustomDate::to_string() const {
    std::string str_rep = std::to_string(year);

    if (month < 10) {
        str_rep += ("0" + std::to_string(month));
    } else {
        str_rep += std::to_string(month);
    }

    if (day < 10) {
        str_rep += ("0" + std::to_string(day));
    } else {
        str_rep += std::to_string(day);
    }

    return str_rep;
}

CustomDate CustomDate::to_bom() const {
    CustomDate date_tmp = *this;
    date_tmp.day = 1;
    date_tmp.hour = 0;
    date_tmp.minute = 0;
    date_tmp.second = 0;

    return date_tmp;
}

CustomDate CustomDate::to_eom() const {
    CustomDate date_tmp = *this;
    date_tmp.day = date_tmp.eom();
    date_tmp.hour = 23;
    date_tmp.minute = 59;
    date_tmp.second = 59;

    return date_tmp;
}

bool CustomDate::is_leap_year(int y) {
    if (((y % 4 == 0) && (y % 100 != 0)) || (y % 400 == 0))
        return true;

    return false;
}

int CustomDate::eom() const {
    int result = max_days[month - 1];
    if (is_leap_year(year) && month == 2) {
        ++result;
    }
    return result;
}

template <typename Comparator>
bool compare(const CustomDate& lhs, const CustomDate& rhs, Comparator comp) {
    const auto lhst = std::make_tuple(
        lhs.get_year(), lhs.get_month(), lhs.get_day(), lhs.get_hour(), lhs.get_minute(), lhs.get_second());
    const auto rhst = std::make_tuple(
        rhs.get_year(), rhs.get_month(), rhs.get_day(), rhs.get_hour(), rhs.get_minute(), rhs.get_second());

    return comp(lhst, rhst);
}

bool operator==(const CustomDate& lhs, const CustomDate& rhs) {
    return compare(lhs, rhs, std::equal_to<std::tuple<int, int, int, int, int, int>>{});
}
bool operator!=(const CustomDate& lhs, const CustomDate& rhs) {
    return compare(lhs, rhs, std::not_equal_to<std::tuple<int, int, int, int, int, int>>{});
}

bool operator<(const CustomDate& lhs, const CustomDate& rhs) {
    return compare(lhs, rhs, std::less<std::tuple<int, int, int, int, int, int>>{});
}

bool operator<=(const CustomDate& lhs, const CustomDate& rhs) {
    return compare(lhs, rhs, std::less_equal<std::tuple<int, int, int, int, int, int>>{});
}

bool operator>(const CustomDate& lhs, const CustomDate& rhs) {
    return compare(lhs, rhs, std::greater<std::tuple<int, int, int, int, int, int>>{});
}

bool operator>=(const CustomDate& lhs, const CustomDate& rhs) {
    return compare(lhs, rhs, std::greater_equal<std::tuple<int, int, int, int, int, int>>{});
}
