#pragma once

#include "LinkedList.h"

class HW_Base : public ListNode
{

public:
    HW_Base();
    virtual ~HW_Base();

    virtual std::string GetName();
    virtual std::string GetNodeType();
    virtual bool SetPropertyMap(std::map<std::string, std::string>* in_map);
    virtual std::map<std::string, std::string>* GetPropertyMap();
    virtual std::string GetProperty(std::string find_name);
    virtual ListNode* clone();

    //+++++Define standard methods for HW_Base - Start+++++
    virtual int GetWidth();
    virtual int GetHeight();
    virtual bool InitMainBoard(LinkedList* m_LinkedList);
    virtual bool UninitMainBoard(LinkedList* m_LinkedList);
    virtual bool RegRead(uchar in_slave_id, uint in_addr, uchar in_addr_len, ushort* in_data, uchar in_data_len);
    virtual bool RegWrite(uchar in_slave_id, uint in_addr, uchar in_addr_len, ushort in_data, uchar in_data_len);
    //++++++Define standard methods for HW_Base - END++++++

    //++++++Define middle methods for HW_Base - Start++++++
    //virtual bool CaptureOriImage(cv::Mat* out_img);
    virtual bool CaptureOriImage(BYTE* buffer, ULONG& grabsize, ushort& uWidth, ushort& uHeight, int& type);    //For Raw Image
    //+++++++Define middle methods for HW_Base - END+++++++
};
