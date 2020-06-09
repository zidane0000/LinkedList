#include "EX_Base.h"

EX_Base::EX_Base()
{
}

EX_Base::~EX_Base()
{

}

std::string EX_Base::GetName()
{
	return "EX_Base";
}

std::string EX_Base::GetNodeType()
{
	return "EX_Base";
}

std::string EX_Base::GetProperty(std::string find_name)
{
	return "";
}

bool EX_Base::PreProcess(LinkedList* m_LinkedList, HW_Base* in_hardware)
{
	//This function isn't complete
	return false;
}

bool EX_Base::Process(LinkedList* m_LinkedList, HW_Base* in_hardware)
{
	//This function isn't complete
	return false;
}

bool EX_Base::PostProcess(LinkedList* m_LinkedList, HW_Base* in_hardware)
{
	//This function isn't complete
	return false;
}

