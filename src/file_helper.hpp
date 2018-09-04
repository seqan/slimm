using namespace seqan;

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <utility>

#ifdef _WIN32
    #include <io.h>
    #define access    _access_s

    std::vector<std::string> get_bam_files_in_directory(std::string directory)
    {
        std::vector<std::string>  input_paths;
        HANDLE dir;
        WIN32_FIND_DATA file_data;

        if ((dir = FindFirstFile((directory + "/*").c_str(),
           &file_data)) == INVALID_HANDLE_VALUE)
            return input_paths; /* No files found */

        do
        {
            const std::string file_name = file_data.cFileName;
            const std::string full_file_name = directory + "/" + file_name;
            const bool is_directory =
            (file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;

            if (file_name[0] == '.')
                continue;

            if (is_directory)
                continue;

            if(full_file_name.find(".sam") == full_file_name.find_last_of(".") ||
                full_file_name.find(".bam")  == full_file_name.find_last_of(".") )
             input_paths.push_back(full_file_name);
        } while (FindNextFile(dir, &file_data));

        FindClose(dir);
        return input_paths;
    }

#else
    #include <unistd.h>
    std::vector<std::string> get_bam_files_in_directory(std::string directory)
    {
        std::vector<std::string>  input_paths;
        DIR *dir;
        struct dirent *ent;
        struct stat st;

        dir = opendir(directory.c_str());
        while ((ent = readdir(dir)) != NULL)
        {
            const std::string file_name = ent->d_name;
            const std::string full_file_name = directory + "/" + file_name;

            if (file_name[0] == '.')
                continue;

            if (stat(full_file_name.c_str(), &st) == -1)
                continue;

            const bool is_directory = (st.st_mode & S_IFDIR) != 0;

            if (is_directory)
                continue;


            if(full_file_name.find(".sam") == full_file_name.find_last_of(".") ||
               full_file_name.find(".bam")  == full_file_name.find_last_of(".") )
                input_paths.push_back(full_file_name);
        }
        closedir(dir);
        return input_paths;
    } // get_bam_files_in_directory

#endif

bool is_file(const char* path)
{
    return access(path, 0 ) == 0;
}

std::string get_file_name (const std::string& str)
{
    std::size_t found = str.find_last_of("/\\");
    return str.substr(found+1);
}

std::string get_directory (const std::string& str)
{
    std::size_t found = str.find_last_of("/\\");
    return str.substr(0,found);
}

std::string get_tsv_file_name(const std::string & output_prefix, const std::string& input_path)
{
    std::string dir_name = get_directory(output_prefix);
    std::string file_name = get_file_name(output_prefix);
    if (file_name.size() == 0)
    {
        file_name = get_file_name(input_path);

        if ((file_name.find(".sam") != std::string::npos &&
           file_name.find(".sam") == file_name.find_last_of("."))
            ||
            (file_name.find(".bam") != std::string::npos &&
               file_name.find(".bam")  == file_name.find_last_of(".")))
        {
            file_name.replace((file_name.find_last_of(".")), 4, "");
        }
    }
    return dir_name + "/" + file_name;
}

std::string get_tsv_file_name (const std::string& output_prefix, const std::string& input_path, const std::string& decor_suffix)
{
    return get_tsv_file_name(output_prefix, input_path) + decor_suffix + ".tsv";
}
