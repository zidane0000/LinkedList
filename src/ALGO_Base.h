#pragma once

#include "LinkedList.h"
#include "HW_Base.h"

class ALGO_Base : public ListNode
{

public:
    ALGO_Base();
    virtual ~ALGO_Base();

    virtual std::string GetName();
    virtual std::string GetNodeType();
    virtual bool GetPassFail();
    virtual bool SetPropertyMap(std::map<std::string, std::string>* in_map);
    virtual std::map<std::string, std::string>* GetPropertyMap();
    virtual std::string GetProperty(std::string find_name);
    virtual ListNode* clone();

    //+++++Define standard methods for ALGO_Base - Start+++++
    virtual bool PreProcess(LinkedList* m_LinkedList, HW_Base* in_hardware);
    virtual bool Process(LinkedList* m_LinkedList, HW_Base* in_hardware);
    virtual bool PostProcess(LinkedList* m_LinkedList, HW_Base* in_hardware);
    //++++++Define standard methods for ALGO_Base - END++++++
};
