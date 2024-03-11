#include "output_helper.h"
#include "scope_timer.hpp"

void field_output(std::vector<std::vector<double>>& field, std::string root_dir, std::string label, int ind)
{
    BeginScopeTimer("field_output(...)");

    std::stringstream ss;
    ss << root_dir << "/" << label << "/";

    if (!fs::exists(ss.str()))
    {
        std::error_code err;
        bool            create_root_dir_okey = fs::create_directory(ss.str(), err);
        std::cout << "create_directory(" << std::quoted(ss.str()) << "), result = " << create_root_dir_okey << "\n";
        std::cout << "err.value() = " << err.value() << "\n"
                  << " err.message() = " << err.message() << "\n";
    }

    ss << label << "_" << ind << ".csv";

    std::string output_filename = ss.str();

    std::ofstream outfile;
    outfile.open(output_filename, std::ios::out);
    for (int i = 0; i < field.size(); i++)
    {
        for (int j = 0; j < field[0].size(); j++)
        {
            if (j == (field[0].size() - 1))
            {
                outfile << field[i][j] << std::endl;
            }
            else
            {
                outfile << field[i][j] << ",";
            }
        }
    }
    outfile.close();
}