// 의생프 - 가변 데이터 처리 문제 연습
#include<iostream>

using namespace std;

int main(){
    int n;

    cout << "입력할 숫자의 개수?" << endl;
    cin >> n;

    // 동적으로! 메모리 할당
    int* numbers = new int[n];

    cout<< n << "개의 숫자를 입력하세요." << endl;

    for (int i = 0; i < n; i++) {
        cin >> numbers[i];
    }

    cout << "입력된 숫자는 다음과 같습니다." << endl;
    for (int i = 0; i < n; i++) {
        cout << numbers[i] << " ";
    }
    cout << endl;

    // 메모리 해제
    delete[] numbers;
}