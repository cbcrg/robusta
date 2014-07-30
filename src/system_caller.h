#ifndef SYSTEM_CALLER_H
#define SYSTEM_CALLER_H

#include "errno.h"        //Error information
#include "fcntl.h"        //FILE CONTROL OPTIONS
#include "stdio.h"
// #include "sys/types.h" //Types
// #include "sys/stat.h"     //File statistics
#include "sys/sysinfo.h"  //Processor control
#include "sys/wait.h"     //Waiting declaration
#include "unistd.h"






/*! \file system_caller.h
	\brief Header file for functions dealing with the handling of system calls.
*/




/**
 * \brief Similar to standard C system().
 *
 * \param arguments The command to execute.
 * \return 0 if no error occured.
 */
int
my_system(char * const arguments[]);




/**
* \brief Similar to C system() but does not wait till the process is finished.
*
* The function does not return immediatly when too many processes are started.
* \param arguments The command to execute.
* \return 0 if no error occured or no child has finished.
*/
int
my_system_no_wait(char **arguments);




/**
 * \brief Similar to C system() but does not wait till the process is finished.
 *
 * This function behaves like \a my_system_no_wait but allows the redirection of stdout and stderr. If stdout or stderr is set to
 * /dev/null the corresponding output will be send into "nothing".
 * \param arguments The command to execute.
 * \param std_out The file in which the standard output should be redirected.
 * \param std_err The file in which the standard error should be redirected.
 * \see my_system_no_wait
 * \see my_system
 */
int
my_system_no_wait_redirect(char **arguments, char *std_out, char *std_err);


/**
* \brief Similar to C system(), waits till the process is finished.
*
* This function behaves like \a my_system_no_wait but allows the redirection of stdout and stderr. If stdout or stderr is set to
* /dev/null the corresponding output will be send into "nothing".
* \param arguments The command to execute.
* \param std_out The file in which the standard output should be redirected.
* \param std_err The file in which the standard error should be redirected.
* \see my_system_no_wait
* \see my_system
*/
int
my_system_redirect(char **arguments, char *std_out, char *std_err);



/**
 * \brief Sets the maximal number of processes to lunch.
 */
void
set_max_num_p(int max_num);

/**
 * \brief Stops the program until all child processes are finished.
 * \return 0 if in none of the child processes an error occured.
 */
int
wait_for_all();


#endif /* SYSTEM_CALLER_H */