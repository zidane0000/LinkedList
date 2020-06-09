#include "HW_Base.h"

HW_Base::HW_Base()
{
}

HW_Base::~HW_Base()
{

}

std::string HW_Base::GetName()
{
	return "HW_Base";
}

std::string HW_Base::GetNodeType()
{
	return "HW_Base";
}

bool HW_Base::SetPropertyMap(std::map<std::string, std::string>* in_map)
{
	return false;
}

std::string HW_Base::GetProperty(std::string find_name)
{
	return "";
}

ListNode* HW_Base::clone()
{
	return nullptr;
}

int HW_Base::GetWidth()
{
	return 0;
}

int HW_Base::GetHeight()
{
	return 0;
}

bool HW_Base::InitMainBoard(LinkedList* m_LinkedList)
{
	//This function isn't complete
	return false;
}

bool HW_Base::UninitMainBoard(LinkedList* m_LinkedList)
{
	//This function isn't complete
	return false;
}

bool HW_Base::RegRead(uchar in_slave_id, uint in_addr, uchar in_addr_len, ushort* in_data, uchar in_data_len)
{
	return false;
}

bool HW_Base::RegWrite(uchar in_slave_id, uint in_addr, uchar in_addr_len, ushort in_data, uchar in_data_len)
{
	return false;
}

//bool HW_Base::CaptureOriImage(cv::Mat* out_img)
//{
//	//This function isn't complete
//	return false;
//}

bool HW_Base::CaptureOriImage(BYTE* buffer, ULONG& grabsize, ushort& uWidth, ushort& uHeight, int& type)
{
	//This function isn't complete
	return false;
}
