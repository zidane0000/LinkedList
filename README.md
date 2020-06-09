# LinkedList 專案介紹
Basic LinkedList , ListNode and Factory

## 介紹
- LinkedList    鏈結串列
- ListNode      節點
- SimpleFactory 工廠模式

## LinkedList
- Push_front(ListNode* newNode);
  - 在list的開頭新增node
- Push_back(ListNode* newNode);
  - 在list的尾巴新增node
- Clear();
  - 把整串list刪除
- getFirst() { return first; };
- getLast() { return last; };
- getNodeByName(std::string NodeName);
- getNodeByValue(int NodePosi);
- Clone();
  - 需先實現ListNode的clone
- HelpRunListNode(ListNode* in_ListNode);

## ListNode
- GetName();
- GetNodeType() { return "ListNode"; };
- SetPropertyMap(std::map< std::string, std::string>* in_map) = 0;
- GetPropertyMap() = 0;
- GetProperty(std::string find_name) = 0;
- clone() = 0;

## Factory
- createListNode(std::string in_name, std::map<std::string, std::string>* new_map)
