
#include "system_caller.h"

static int max_num_p = 0;
static int current_num_p = 0;

void
set_max_num_p(int max_num)
{
	if (max_num > 0)
		max_num_p = max_num;
	else
		max_num_p = get_nprocs ();
}



int
my_system(char *const arguments[])
{
	int status;
	++current_num_p;
	pid_t pid = vfork();

	if (pid < 0)
	{
		_exit(2);
	}
	else if (pid == 0)
	{
		setpgid(0, getppid());
		execvp(arguments[0], arguments);

		//Since exevcp only returns when not successfull this code is only executed when an error occured.
		_exit(1);
	}
	else
	{
		waitpid(pid, &status, WUNTRACED);
		--current_num_p;
		return status;

	}
}



int
my_system_redirect(char **arguments, char *std_out, char *std_err)
{
	++current_num_p;
	pid_t pid = vfork();

	if (pid < 0)
	{
		_exit(2);
	}
	else if (pid == 0)
	{
		setpgid(0, getppid());

		int fd;
		if((fd = open(std_out, O_RDWR  , S_IWUSR | S_IRUSR))==-1)
		{
			perror("open");
			return 1;
		}
		dup2(fd, STDOUT_FILENO);
		close(fd);
		if((fd = open(std_err, O_RDWR | O_CREAT , S_IWUSR | S_IRUSR))==-1)
		{
			perror("open");
			return 1;
		}
		dup2(fd, STDERR_FILENO);
		close(fd);
		// 		printf(arguments[0], arguments[1], arguments
		execvp(arguments[0], arguments);

		//Since exevcp only returns when not successfull this code is only executed when an error occured.
		_exit(1);
	}
	else
	{
		int status;
		waitpid(pid, &status, WUNTRACED);
		--current_num_p;
		return status;
	}

}





int
my_system_no_wait(char **arguments)
{
	++current_num_p;
	pid_t pid = vfork();

	if (pid < 0)
	{
		_exit(2);
	}
	else if (pid == 0)
	{
		setpgid(0, getppid());
		execvp(arguments[0], arguments);

		//Since exevcp only returns when not successfull this code is only executed when an error occured.
		_exit(1);
	}
	else
	{

		if (current_num_p >= max_num_p)
		{
			int status;
			wait(&status);
			--current_num_p;
			return status;
		}
		return 0;
	}

}


int
my_system_no_wait_redirect(char **arguments, char *std_out, char *std_err)
{
	++current_num_p;
	pid_t pid = vfork();

	if (pid < 0)
	{
		_exit(2);
	}
	else if (pid == 0)
	{
		setpgid(0, getppid());

		int fd;
		if((fd = open(std_out, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR))==-1)
		{
			perror("open");
			return 1;
		}
		dup2(fd, STDOUT_FILENO);
		close(fd);
		if((fd = open(std_err, O_RDWR | O_CREAT , S_IWUSR | S_IRUSR))==-1)
		{
			perror("open");
			return 1;
		}
		dup2(fd, STDERR_FILENO);
		close(fd);

		execvp(arguments[0], arguments);

		//Since exevcp only returns when not successfull this code is only executed when an error occured.
		_exit(1);
	}
	else
	{

		if (current_num_p >= max_num_p)
		{
			int status;
			wait(&status);
			--current_num_p;
			return status;
		}
		return 0;
	}

}




int
wait_for_all()
{
	pid_t pid;
	int status, total = 0;
	while ((pid = waitpid(-1, &status, 0))) {
		total += status;
		if (errno == ECHILD) { //EXIT when no child exists anymore
			break;
		}
		--current_num_p;
	}
	return total;
}

