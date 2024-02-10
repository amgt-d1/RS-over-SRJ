#include "isrjs.hpp"
#include <time.h>


int main()
{
    // input parameters
    input_parameter();

    // input data
    input_data();

    time_t t = time(NULL);
	printf(" %s\n\n", ctime(&t));

    isrjs ISRJS;

    // join sampling
    ISRJS.join_sampling();

    return 0;
}