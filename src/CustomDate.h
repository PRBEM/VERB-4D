/*
 * SPDX-FileCopyrightText: 2015 UCLA
 * SPDX-FileCopyrightText: 2025 Bernhard Haas (GFZ)
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/**
 * \file CustomDate.h
 * \brief Custom date and time handling class
 */

#pragma once
#include <ctime>
#include <string>

/**
 * @brief Custom date class with MATLAB datenum compatibility
 */
class CustomDate {
   public:
    /** @brief MATLAB datenum epoch value */
    static const double matlabd_epoch;
    
    /** @brief Default constructor */
    CustomDate() {}
    
    /** @brief Constructs date from MATLAB datenum value */
    CustomDate(double matlabd);
    
    /** @brief Copy constructor */
    CustomDate(const CustomDate& rhs);
    
    /** @brief Assignment operator */
    CustomDate& operator=(const CustomDate& rhs);
    
    /** @brief Adds seconds to the date */
    CustomDate operator+(uint32_t delta_seconds);
    
    /** @brief Subtracts seconds from the date */
    CustomDate operator-(uint32_t delta_seconds);
    
    /** @brief Advances date by one month */
    CustomDate add_month() const;

    /** @brief Converts date to string format YYYYMMDD */
    std::string to_date_string() const;
    
    /** @brief Converts date and time to ISO format string */
    std::string to_string() const;
    
    /** @brief Returns date at beginning of current month */
    CustomDate to_bom() const;
    
    /** @brief Returns date at end of current month */
    CustomDate to_eom() const;
    
    /** @brief Equality comparison operator */
    bool operator==(const CustomDate& rhs);
    
    /** @brief Inequality comparison operator */
    bool operator!=(const CustomDate& rhs);
    
    /** @brief Less than comparison operator */
    bool operator<(const CustomDate& rhs);
    
    /** @brief Less than or equal comparison operator */
    bool operator<=(const CustomDate& rhs);
    
    /** @brief Greater than comparison operator */
    bool operator>(const CustomDate& rhs);
    
    /** @brief Greater than or equal comparison operator */
    bool operator>=(const CustomDate& rhs);
   private:
    /** @brief Standard time structure for date/time components */
    std::tm datetime;
    
    /** @brief Unix timestamp in seconds */
    std::time_t unix_seconds = 0;
    
    /** @brief Maximum days in each month */
    static const int max_days[12];
    
    /** @brief Checks if current year is a leap year */
    bool is_leap_year() const;
    
    /** @brief Returns last day of current month */
    int eom() const;
};
