#ifndef LLALBM_LOGGER_HPP
#define LLALBM_LOGGER_HPP

// =========== STL INCLUDES ===========
#include <iostream>
#include <ostream>
#include <string>
// ======================================

namespace llalbm::util::logger
{
    class Logger
    {
    private:
        /// @brief Name of the logger, used while printing messages
        const std::string name;
        /// @brief Output stream onto which the message will be printed
        std::ostream& out_stream;
    public:
        Logger(const std::string& name_, std::ostream& out_stream_)
        :   name (name_)
        ,   out_stream (out_stream_)
        {}

        inline void info_no_return(const std::string& msg)
        {
            out_stream << "["<<name<<"-INFO] "<< msg;
        }

        inline void info(const std::string& msg)
        {
            out_stream << "["<<name<<"-INFO] "<< msg << std::endl;
        }

        inline void warn(const std::string& msg)
        {
            out_stream << "["<<name<<"-WARNING] "<< msg << std::endl;
        }

        inline void error(const std::string& msg)
        {
            out_stream << "["<<name<<"-ERROR] "<< msg << std::endl;
        }
    };
};

#endif