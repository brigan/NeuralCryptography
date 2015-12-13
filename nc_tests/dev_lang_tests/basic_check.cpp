//
// Created by Nurlan Akashayev on 12/12/15.
//

#include <iostream>
#include <list>
#include "gtest/gtest.h"
using namespace std;

TEST(basic_check, test_list_deletion) {

    std::list<int> list1;

    // Add some values at the end of the list, which is initially empty.
    // The member function "push_back" adds at item at the end of the list.
    int value1 = 10;
    int value2 = -3;
    list1.push_back (value1);
    list1.push_back (value2);
    list1.push_back (5);
    list1.push_back (1);

    std::list<int>::iterator i = list1.begin();
//    while (i != list1.end()) {
//        //cout << *i;
//    }

    EXPECT_EQ(4, list1.size());
}

