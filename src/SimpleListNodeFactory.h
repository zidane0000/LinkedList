#pragma once

#include "LinkedList.h"
#include "HW_Base.h"
#include "ALGO_Base.h"
#include "EX_Base.h"

class SimpleListNodeFactory
{

public:
	ListNode* createListNode(std::string in_name, std::map<std::string, std::string>* new_map) {
		ListNode* new_node = nullptr;

		//+++Create ListNode - Start+++
		if (in_name.find("HW_") != in_name.npos)
		{
			new_node = new HW_Base();
		}
		else if (in_name.find("ALGO_") != in_name.npos)
		{
			new_node = new ALGO_Base();
		}
		else if (in_name.find("EX_") != in_name.npos)
		{
			new_node = new EX_Base();
		}
		//++++Create ListNode - END++++

		if(new_node)
			new_node->SetPropertyMap(new_map);
		return new_node;
	}

};
