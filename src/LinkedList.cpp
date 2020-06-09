#include "LinkedList.h"
#include "HW_Base.h"
#include "ALGO_Base.h"
#include "EX_Base.h"

void LinkedList::Push_front(ListNode* newNode) {
    newNode->next = first;                  // 先把first接在newNode後面
    first = newNode;                        // 再把first指向newNode所指向的記憶體位置
}

void LinkedList::Push_back(ListNode* newNode){
    last->next = newNode;                  // 先把newNode接在last後面
    last = newNode;                        // 再把last指向newNode所指向的記憶體位置
}

void LinkedList::Clear()
{
    while (first != 0) {            // Traversal
        ListNode* current = first;
        first = first->next;
        delete current;
        current = 0;
    }
    first = nullptr;
    last = nullptr;
}

ListNode* LinkedList::getNodeByName(std::string NodeName)
{
    ListNode* current_node = first;
    while (current_node != NULL && (current_node->GetName() != NodeName)) {
        current_node = current_node->next;
    }

    if (current_node->GetName() == NodeName)
        return current_node;
    else
        return nullptr;
}

ListNode* LinkedList::getNodeByValue(int NodePosi)
{
    ListNode* current_node = first;
    while (NodePosi-- && (current_node != nullptr)) {
        current_node = current_node->next;
    }
    return current_node;
}

LinkedList* LinkedList::Clone()
{
    LinkedList* new_LinkedList = new LinkedList();

    ListNode* current_node = first;
    ListNode* current_clone = current_node->clone();
    new_LinkedList->SetFirstAndLast(current_clone);

    current_node = current_node->next;
    while (current_node != NULL) {
        current_clone = current_node->clone();
        new_LinkedList->Push_back(current_clone);
        current_node = current_node->next;
    }

	return new_LinkedList;
}

std::string ListNode::GetName()
{
    //This function isn't complete
    return "ListNode";
}

bool LinkedList::HelpRunListNode(ListNode* in_ListNode)
{
    std::string NodeTypeName = in_ListNode->GetNodeType();
    int nRes = -1;
    if (NodeTypeName.find("HW_") != NodeTypeName.npos) {
        std::string NodeName = in_ListNode->GetName();
        
        HW_map[NodeName] = in_ListNode;
    }
    else if (NodeTypeName.find("ALGO_") != NodeTypeName.npos) {
        ALGO_Base* at_ALGO = dynamic_cast<ALGO_Base*>(in_ListNode);

        //+++Find HW - Start +++
        std::string at_HW_name = at_ALGO->GetProperty("HW_Name");
        HW_Base* at_HW = dynamic_cast<HW_Base*>(HW_map[at_HW_name]);
        if (!at_HW) {
            ThreadisPass = false;
            return nRes;
        }
        //++++Find HW - END ++++

        //+++ Run at_ALGO - Start +++
        nRes = at_ALGO->PreProcess(this, at_HW);
        nRes = at_ALGO->Process(this, at_HW);
        nRes = at_ALGO->PostProcess(this, at_HW);
        //++++ Run at_ALGO - END ++++

        //Result Show On UI
        std::string NodeName = at_ALGO->GetName();
        if (at_ALGO->GetPassFail()) {
            
        }
        else {
            ThreadisPass = false;
        }
    }
    else if (NodeTypeName.find("EX_") != NodeTypeName.npos) {
        EX_Base* at_EX = dynamic_cast<EX_Base*>(in_ListNode);

        //+++Find HW - Start +++
        std::string at_HW_name = at_EX->GetProperty("HW_Name");
        HW_Base* at_HW = dynamic_cast<HW_Base*>(HW_map[at_HW_name]);
        if (!at_HW) {
            ThreadisPass = false;
            return nRes;
        }
        //++++Find HW - END ++++

        //+++ Run at_EX - Start +++
        nRes = at_EX->PreProcess(this, at_HW);
        nRes = at_EX->Process(this, at_HW);
        nRes = at_EX->PostProcess(this, at_HW);
        //++++ Run at_EX - END ++++

         //Result Show On UI
        std::string NodeName = at_EX->GetName();
        if (at_EX->GetPassFail()) {

        }
        else {
            ThreadisPass = false;
        }
    }
    else {

    }

    return nRes;
}