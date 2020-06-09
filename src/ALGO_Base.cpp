#include "ALGO_Base.h"

ALGO_Base::ALGO_Base()
{
}

ALGO_Base::~ALGO_Base()
{

}

std::string ALGO_Base::GetName()
{
	return "ALGO_Base";
}

std::string ALGO_Base::GetNodeType()
{
	return "ALGO_Base";
}

std::string ALGO_Base::GetProperty(std::string find_name)
{
	return "";
}

bool ALGO_Base::PreProcess(LinkedList* m_LinkedList, HW_Base* in_hardware)
{
	//This function isn't complete
	return false;
}

bool ALGO_Base::Process(LinkedList* m_LinkedList, HW_Base* in_hardware)
{
	//This function isn't complete
	return false;
}

bool ALGO_Base::PostProcess(LinkedList* m_LinkedList, HW_Base* in_hardware)
{
	//This function isn't complete
	return false;
}

