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

std::string HW_Base::GetProperty(std::string find_name)
{
	return "";
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
