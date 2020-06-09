#pragma once
#include <iostream>
#include <string>
#include <map>
#include <vector>

typedef unsigned        uint;
typedef unsigned char   uchar;
typedef unsigned short  ushort;
typedef unsigned char   BYTE;
typedef unsigned long   ULONG;

class LinkedList;    // ���F�Nclass LinkedList�]��class ListNode��friend,�ݭn���ŧi

class ListNode
{

public:
    virtual std::string GetName();
    virtual std::string GetNodeType() { return "ListNode"; };
    virtual bool SetPropertyMap(std::map< std::string, std::string>* in_map) = 0;
    virtual std::map< std::string, std::string>* GetPropertyMap() = 0;
    virtual std::string GetProperty(std::string find_name) = 0;
    virtual ListNode* clone() = 0;
    
    friend class LinkedList;
    ListNode* next = nullptr;
};

class LinkedList
{
    
public:
    LinkedList() {
        first = nullptr; last = nullptr;
    };

    bool SetFirstAndLast(ListNode* newNode) {
        if ((first == nullptr) && (last == nullptr)) {
            first = newNode; last = newNode;
            return true;
        }
        else
            return false;
    };

    void Push_front(ListNode* newNode);     // �blist���}�Y�s�Wnode
    void Push_back(ListNode* newNode);      // �blist�����ڷs�Wnode
    void Clear();               // ����list�R��
    ListNode* getFirst() { return first; };
    ListNode* getLast() { return last; };
    ListNode* getNodeByName(std::string NodeName);
    ListNode* getNodeByValue(int NodePosi);
    LinkedList* Clone();
    bool HelpRunListNode(ListNode* in_ListNode);

private:
    ListNode* first;            // list���Ĥ@��node
    ListNode* last;             // list���̫�@��node
    
    int node_posi = 0;          //�{�b�B������
    bool ThreadisPass = false;

    std::map<std::string, ListNode*> HW_map;
};
