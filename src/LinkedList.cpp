#include "LinkedList.h"
void LinkedList::Push_front(ListNode* newNode) {
    newNode->next = first;                  // ����first���bnewNode�᭱
    first = newNode;                        // �A��first���VnewNode�ҫ��V���O�����m
}

void LinkedList::Push_back(ListNode* newNode){
    last->next = newNode;                  // ����newNode���blast�᭱
    last = newNode;                        // �A��last���VnewNode�ҫ��V���O�����m
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