#ifndef DEF_INPUT
#define DEF_INPUT

#include <string>
#include <vector>

void load_CSV(std::string base_name, bool use_CNA=true, bool apply_filter_regions=true);
void filter_regions();
void init_params();

#endif