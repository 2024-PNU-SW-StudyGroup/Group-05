// 의생프 - 전역변수 지역변수 연습
#include <iostream>
using namespace std;

int globalVar = 10; // 전역 변수
void exampleFunction() {
 int localVar = 20; // 함수의 지역 변수

 cout << "전역 변수: " << globalVar << endl;
 cout << "함수의 지역 변수: " << localVar << endl;

 {
 int nestedVar = 30;
 cout << "중첩된 블록의 지역 변수: " << nestedVar << endl;

 int globalVar = 40;
 cout << "같은 이름의 지역 변수: " << globalVar << endl;
 cout << "전역 변수 (::globalVar): " << ::globalVar << endl;
 }

 // cout << nestedVar << endl; // 오류: nestedVar는 여기서 접근할 수 없음
}

int main() {
 exampleFunction();

 // cout << localVar << endl; // 오류: localVar는 여기서 접근할 수 없음

 return 0;
}