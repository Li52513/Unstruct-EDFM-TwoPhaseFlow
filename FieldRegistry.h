#pragma once
#include <unordered_map>
#include <memory>
#include <string>
#include <utility>
#include "VolField.h"

// FieldRegistry 可以实现以下功能：  
// 1. 动态创建物理场对象并存储：通过 `create` 方法，可以动态创建任意类型的物理场对象（继承自 BaseField），并将其存储在字段映射中。  
// 2. 获取已存储的物理场对象：通过 `get` 方法，可以根据名称获取存储的物理场对象，并支持动态类型转换。  
// 3. 检查物理场对象是否存在：通过 `has` 方法，可以检查某个名称的物理场对象是否已经存在于映射中。  

struct FieldRegistry 
{
	unordered_map<string, shared_ptr<BaseField>> fields; //建立名称到物理场变量的映射

template<class FieldType, class... Args>
shared_ptr<FieldType> create(const string& name, Args&&... args) 
{
    // 使用 std::make_shared 创建一个 FieldType 对象
    auto f = std::make_shared<FieldType>(name, forward<Args>(args)...);

    // 将创建的对象存储到 fields 映射中
    fields[name] = f;

    // 返回智能指针，供调用者使用
    return f;
}


    template<class FieldType>
    shared_ptr<FieldType> get(const string& name) const {
        auto it = fields.find(name);
        if (it == fields.end()) return nullptr;
        return dynamic_pointer_cast<FieldType>(it->second);
    }

    bool has(const string& name) const { return fields.count(name) > 0; }
};