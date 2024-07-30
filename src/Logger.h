/**
 * \file Logger.h
 * \brief Logs the messages sent, including information, warnings, and errors
 */

#pragma once

#include <iostream>
#include <streambuf>
#include <fstream>
#include <sstream>
#include <stdio.h> 
#include <assert.h>

static const std::string LOGGER_SEPARATOR = std::string(100, '=');

/**
* \brief Stores messages about statuses, warning, and errors when computing/solving the matrices
*/
class Logger
{
public:

	static void writeSeparator() {
		Logger::message << LOGGER_SEPARATOR << std::endl;
	}
	
/**
* Enum for the Debug Level that the logger is at
*/
	enum DebugLevel
	{
		DEBUG_LEVEL_DISABLED = 0, /**< Disabled = 0 */
		DEBUG_LEVEL_MESSAGE = 1, /**< Message = 1 */
		DEBUG_LEVEL_WARNING = 2, /**< Warning = 2 */
		DEBUG_LEVEL_ERROR = 3,	/**< Error = 3 */
	};

/**
* Enum for the type of Message being recorded
*/
	enum MessageType
	{
		MESSAGE_INFO = 1, /**< Info = 1 */
		MESSAGE_WARNING = 2, /**< Warning = 2 */
		MESSAGE_ERROR = 3, /**< Error = 3 */
	};

	static void createInstance();
	static void deleteInstance();
	static void setDebugLevel(DebugLevel level);

/**
* \brief Turns messages from the logger into streams
*/
	class Streamer : public std::ostream {
	public:
		Streamer(Logger::MessageType messageType);
		~Streamer();

	private:
		/**
		* \brief Creates a buffer for use by the Streamer
		*/
		class StringBuffer : public std::stringbuf {
		public:
			StringBuffer(Logger::MessageType messageType);
			~StringBuffer();
			virtual int sync();

		private:
			Logger::MessageType mMessageType;
		};
	};

	static Streamer message;
	static Streamer warning;
	static Streamer error;

	

protected:
	/// defualt constructor
	Logger() {};
	/// copy constructor
	Logger(const Logger& thisobject);
	/// overloaded equal for copy constructor
	Logger& operator=(const Logger&);
	/// default destructor
	~Logger() {};

private:
	/// private member function here
	void writeMessage(const std::string& message);
	void writeWarning(const std::string& message);
	void writeError(const std::string& message);

	void appendToFile(const std::string& message);

	static DebugLevel mDebugLevel;
	static Logger* sInstance; // ??
	static std::ofstream logFile;
};


	

	//Logger& operator<< (std::ostream& (*pfun)(std::ostream&));
	/*Logger& operator<<(const char* sb);
	friend Logger operator<<(Logger &log, int i);*/
	
	//template <typename ANY>
	//friend Logger& operator<<(Logger& logClass, ANY val);
	
	
	/*friend Logger& operator<<(Logger &log, int i);
	friend Logger& operator<<(Logger &log, long i);
	friend Logger& operator<<(Logger &log, float i);
	friend Logger& operator<<(Logger &log, double i);
	friend Logger& operator<<(Logger &log, const char* sb);
	friend Logger& operator<<(Logger &log, string sb);

	template <class T>
	friend Logger& operator<<(Logger &log, T i);
	Logger& operator << (std::ostream& (*pfun)(std::ostream&));*/
//};

