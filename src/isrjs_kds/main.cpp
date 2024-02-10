#include "isrjs_kds.hpp"


int main()
{
    // input parameters
    input_parameter();

    // input data
    input_data();

    time_t t = time(NULL);
	printf(" %s\n\n", ctime(&t));

    // join sampling
    isrjs_kds ISRJS;
    ISRJS.join_sampling();

    return 0;
}