#pragma once

#ifdef NDEBUG
#define DBG(x)
#define DBG_MSG(x)
#define DBG_SET_LOG_FILE(x)
#define DBG_HERE(ignoredArgument)
#define DBG_PEXPR(x)
#define DBG_PVAL(x)
#define DBG_HALT(ignoredArgument)

#else

#define DBG(x) do { x; } while(false)

#include <iostream>
#include <cstdint>
#include <utility>
#include <map>
#include <string>
#include <iterator>
#include <fstream>
#include <cstdlib>

namespace dbg {

   static std::ofstream logFile;
   static std::ostream* dbgOut = &std::clog;

   inline void
   setLogFile(
      const char* fileName
   ) {
      static bool redirectedToLog(false);
      if (redirectedToLog) {
         logFile.close();
      }

      logFile.open(fileName);

      dbgOut = (std::ostream*) &logFile;
      redirectedToLog = true;
   }

   #define DBG_SET_LOG_FILE(fileName) DBG(::dbg::setLogFile(fileName))


   #define DBG_MSG(message) DBG((*::dbg::dbgOut) << __FILE__ ":" << __LINE__ << " " << message << "\n")


   inline void
   here(
      const std::string& fileName,
      const int lineNo
   ) {
      static std::map<std::pair<std::string, int>, unsigned> timesBeenAt;
      const std::pair<std::string, int> location(std::make_pair(fileName, lineNo));

      unsigned count;
      std::map<std::pair<std::string, int>, unsigned>::iterator currentLocation(timesBeenAt.find(location));
      if (currentLocation == timesBeenAt.end()) {
         count = timesBeenAt[location] = 1;
      } else {
         count = ++(currentLocation->second);
      }

      DBG_MSG("Been here " << count << " times");
   }

   #define DBG_HERE(ignoredArgument) ::dbg::here(__FILE__, __LINE__)

   #define DBG_PEXPR(variable) DBG_MSG(#variable " = " << variable)
   #define DBG_PVAL(variable) DBG_MSG(variable)
   #define DBG_HALT(ignoredArgument) DBG(DBG_MSG("Halting..."); std::abort())
}


#endif
