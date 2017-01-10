//
//  load_mapping_files.h
//  seqan
//
//  Created by Temesgen H. Dadi on 04/09/15.
//
//


#ifndef APP_SLIMM_LOAD_MAPPING_FILES_H_
#define APP_SLIMM_LOAD_MAPPING_FILES_H_

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <utility>

template <typename TTarget, typename TString, typename TKey = uint32_t, typename TValue = uint32_t>
TTarget loadMapping(TString const & filePath)
{
    TTarget result;
    std::ifstream giMap(toCString(filePath));
    TKey key;
    TValue value;
    while (giMap >> key >> value)
    {
        result[key]=value;
    }
    giMap.close();
    return result;
}

template <typename TTarget, typename TString>
TTarget loadMappingInt2String(TString const & filePath)
{
    TTarget result;
    std::ifstream nameMap(filePath);
    uint32_t key;
    std::string value,  line;
    
    while(std::getline(nameMap, line))
    {
        std::stringstream   linestream(line);
        linestream >> key;
        std::getline(linestream, value, '\t');
        std::getline(linestream, value, '\t');
        result[key]=value;
    }
    nameMap.close();
    return result;
}

template <typename TTarget, typename TString, typename TKey = uint32_t,
            typename TValue1 = uint32_t, typename TValue2 = std::string>
void loadNodes(TTarget & target, TString const & filePath)
{
    std::ifstream nodeMap(filePath);
    std::string   line;
    
    while(std::getline(nodeMap, line))
    {
        std::stringstream   linestream(line);
        TKey key;
        TValue1 value1;
        TValue2 value2;
        linestream >> key >> value1;
        std::getline(linestream, value2, '\t');
        std::getline(linestream, value2, '\t');        
        target[key]=std::pair<TValue1, TValue2> (value1, value2);
    }
    
    nodeMap.close();
}
#endif //APP_SLIMM_LOAD_MAPPING_FILES_H_
