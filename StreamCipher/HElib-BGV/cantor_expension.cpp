#include <bits/stdc++.h>
using namespace std;

const int MAX_N = 36;
__int128 jc[MAX_N]; // 存储 0! 到 35! 的阶乘值
int num[MAX_N]; // num[i] 存储从小到大第 i 个未用过的字母编号

void init() {
    for (int i = 0; i < MAX_N; i++)
        num[i] = i + 1;
}

void precompute_factorials() {
    jc[0] = 1;
    for (int i = 1; i < MAX_N; i++)
        jc[i] = jc[i - 1] * i;
}

vector<int> get_permutation(__int128 cnt) { // 逆展开
    vector<int> ans;
    vector<int> available(num, num + MAX_N); // 使用一个局部数组来存储可用的数字
    for (int i = 0; i < MAX_N - 1; i++) {
        int tmp = cnt / jc[MAX_N - 2 - i];
        ans.push_back(available[tmp] - 1);
        available.erase(available.begin() + tmp); // 直接删除已使用的数字
        cnt %= jc[MAX_N - 2 - i];
    }
    ans.push_back(available[0] - 1); // 添加最后一个剩余的数字
    return ans;
}

void print_permutation(const vector<int>& permutation) {
    for (int i : permutation) {
        cout << i << " ";
    }
    cout << endl;
}

int main() {
    precompute_factorials();
    __int128 n = 12345678901234567889ULL; // 固定的整数，0~36!-1
    init();
    clock_t start = clock();
    vector<int> permutation;
    for (int i = 0; i < 32768; i++) {
        permutation = get_permutation(n);
    }
    clock_t end = clock();
    cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
    // print_permutation(permutation);
    return 0;
}