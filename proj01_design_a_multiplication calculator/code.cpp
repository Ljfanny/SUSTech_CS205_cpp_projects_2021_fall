#include <iostream>
#include <cstdint>
#include <vector>
#include <stack>
#include <string>
#include <sstream>

using namespace std;

string karatsubaMul(string n, string m, int zenum);

string intMul(const int a[], const int b[], int lenA, int lenB, string signDot);

int judgeCal(string c, string d);

bool judgeLawful(string c, string d);

string judSignDot(string c, string d);

bool calculation(string a, string b);

string myPlus(string a, string b, string sign);

int main(int argc, char **argv) {
    string a, b;
    //获取命令行参数
    if (argc > 1) {
        a = argv[1];
        b = argv[2];
    } else {
        //从操作行获取参数
        cout << "Please input two numbers: \n";
        cin >> a;
        cin >> b;
    }

    bool calRes = calculation(a, b);

    if (!calRes) {
        while (!judgeLawful(a, b)) {
            cout << "Please input again: \n";
            cin >> a;
            cin >> b;
            if (calculation(a, b)) {
                exit(0);
            }
        }
        int judCalAns = judgeCal(a, b);
        while (judCalAns == 2) {
            cout << "Please input again: \n";
            cin >> a;
            cin >> b;
            if (calculation(a, b))
                exit(0);
            else
                judCalAns = judgeCal(a, b);
        }
        if (judCalAns == 1) {
            cout << "Calculation successful! \n";
            cout << "The answer is: " << a << " * " << b << endl;
        }
        if (judCalAns == 3) {
            string answer = karatsubaMul(a, b, 0);
            int zeroNum = 0;
            int firLen = answer.length();
            for (int i = 0; i < answer.length(); ++i) {
                if (answer[i] == '0')
                    zeroNum++;
                else break;
            }

            cout << "Calculation successful! \n";
            cout << "The answer is:" << answer.substr(zeroNum, firLen) << endl;
        }
    }
}

string intMul(const int a[], const int b[], int lenA, int lenB, string signDot) {
    string sig, dt;
    istringstream is(signDot);
    is >> sig >> dt;
    if (lenA == 1 && a[0] == 0)
        return "0";
    if (lenB == 1 && b[0] == 0)
        return "0";
    string res;
    int matrix[lenB][lenA + lenB];
    int result[lenA + lenB];
    for (int i = 0; i < lenB + lenA; i++) {
        result[i] = 0;
    }
    for (int i = 0; i < lenB; ++i) {
        for (int j = 0; j < lenB + lenA; ++j) {
            matrix[i][j] = 0;
        }
    }
    for (int i = 0; i < lenB; i++) {
        for (int j = 0; j < lenA; j++) {
            int temp = matrix[i][j + i] + b[i] * a[j];
            matrix[i][j + i] = temp % 10;
            matrix[i][j + i + 1] = (temp - temp % 10) / 10;
        }
    }
    for (int i = 0; i < lenB + lenA; i++) {
        int tem = result[i];
        for (int j = 0; j < lenB; j++) {
            tem += matrix[j][i];
        }
        result[i] = tem % 10;
        if (i + 1 < lenA + lenB)
            result[i + 1] = (tem - result[i]) / 10;
    }

    stack<int> temp;
    for (int i = 0; i < lenA + lenB; ++i) {
        temp.push(result[i]);
    }
    while (!temp.empty()) {
        res += to_string(temp.top());
        temp.pop();
    }

    if (dt != "0")
        res = res.insert(res.length() - atoi(dt.c_str()), 1, '.');
    long firInt = 0;
    long firDot = 0;
    long lastDot = 0;
    long lastInt = 0;
    for (int i = 0; i < res.length(); ++i) {
        if (res[i] != '0') {
            if (res[i] == '.')
                firDot = i;
            else
                firInt = i;
            break;
        }
    }
    if (firDot != 0)
        res = res.substr(firDot - 1, res.length());
    else if (firInt != 0)
        res = res.substr(firInt, res.length());

    for (int i = 0; i < res.length(); ++i) {
        if (res[i] == '.') {
            lastDot = i;
            continue;
        }
        if (lastDot != 0 && res[i] != '0')
            lastInt = i;
    }
    if (lastInt > lastDot)
        res = res.substr(0, lastInt + 1);
    else if (lastInt < lastDot)
        res = res.substr(0, lastDot);
    if (sig == "-")
        res = sig + res;
    return res;
}

bool judgeLawful(string c, string d) {

    if ((long) c.length() > INT32_MAX || (long) c.length() < INT32_MIN) {
        cout << "Your input is too long to calculate!\n";
        return false;
    }
    if ((long) d.length() > INT32_MAX || (long) d.length() < INT32_MIN) {
        cout << "Your input is too long to calculate!\n";
        return false;
    }

    for (int i = 0; i < c.length(); ++i) {
        if (c[i] < 48 || (c[i] > 57 && c[i] < 65) || (c[i] > 90 && c[i] < 97) ||
            c[i] > 122) {
            cout << "Your input contains non-numerals or non-letters and it is invalid!\n";
            return false;
        }
    }
    for (int i = 0; i < d.length(); ++i) {
        if (d[i] < 48 || (d[i] > 57 && d[i] < 65) || (d[i] > 90 && d[i] < 97) ||
            d[i] > 122) {
            cout << "Your input contains non-numerals or non-letters and it is invalid!\n";
            return false;
        }
    }
    return true;
}

int judgeCal(string c, string d) {
    int ans = 0;
    string judChar;
    for (int i = 0; i < c.length(); ++i) {
        if (ans == 1 || ans == 2) {
            break;
        }
        if ((c[i] >= 97 && c[i] <= 122) || (c[i] >= 65 && c[i] <= 90)) {
            cout << "Your input contains non-numbers!\n";
            cout << "Do you still want to calculate?\n";
            cout << "Please input y if you want and input n if you not.\n";
            cin >> judChar;
            while (judChar != "y" && judChar != "n") {
                cout << "Your input is wrong! \n";
                cin >> judChar;
            }
            if (judChar == "y") {
                ans = 1;
            } else {
                ans = 2;
            }
        }
    }
    for (int i = 0; i < d.length(); ++i) {
        if ((d[i] >= 97 && d[i] <= 122) || (d[i] >= 65 && d[i] <= 90)) {
            if (ans == 0) {
                cout << "Your input contains non-numbers!\n ";
                cout << "Do you still want to calculate?\n";
                cout << "Please input y if you want and input n if you not.\n";
                cin >> judChar;
                while (judChar != "y" && judChar != "n") {
                    cout << "Your input is wrong!\n";
                    cin >> judChar;
                }
                if (judChar == "y") {
                    ans = 1;
                } else {
                    ans = 2;
                }
            }
        }
        if (ans == 1 || ans == 2)
            break;
    }
    if (ans != 2 && ans != 1)
        ans = 3;
    return ans;
}

string judSignDot(string c, string d) {
    int dotC = 0;
    int dotD = 0;
    int positionC = 0, positionD = 0;
    string ret;

    if (c[0] != '-' && (c[0] < 48 || c[0] > 57))
        return "false";
    if (d[0] != '-' && (d[0] < 48 || d[0] > 57))
        return "false";

    for (int i = 1; i < c.length(); ++i) {
        if (c[i] == '.') {
            dotC++;
        }
        if (dotC == 1 && c[i] == '.') {
            positionC = c.length() - 1 - i;
            continue;
        }
        if (dotC > 1) {
            return "false";
        }
        if (c[i] < 48 || c[i] > 57) {
            return "false";
        }
    }
    for (int i = 1; i < d.length(); ++i) {
        if (d[i] == '.') {
            dotD++;
        }
        if (dotD == 1 && d[i] == '.') {
            positionD = d.length() - 1 - i;
            continue;
        }
        if (dotD > 1) {
            return "false";
        }
        if (d[i] < 48 || d[i] > 57) {
            return "false";
        }
    }
    if ((c[0] == '-' && d[0] == '-') || (c[0] != '-' && d[0] != '-')) {
        ret = "+ " + to_string(positionD + positionC);
    } else if ((c[0] != '-' && d[0] == '-') || (c[0] == '-' && d[0] != '-'))
        ret = "- " + to_string(positionD + positionC);

    return ret;
}

bool calculation(string a, string b) {
    if (judSignDot(a, b) != "false") {
        int countA = 0;
        int countB = 0;
        if (a[0] == '-')
            countA++;
        for (int i = 0; i < a.length(); ++i) {
            if (a[i] == '.')
                countA++;
        }
        if (b[0] == '-')
            countB++;
        for (int i = 0; i < b.length(); ++i) {
            if (b[i] == '.')
                countB++;
        }
        string usA, usB;

        for (int i = a.length() - 1; i >= 0; --i) {
            if (a[a.length() - i - 1] >= 48 && a[a.length() - i - 1] <= 57) {
                usA += to_string(a[a.length() - i - 1] - '0');
            }
        }

        for (int i = b.length() - 1; i >= 0; --i) {
            if (b[b.length() - i - 1] >= 48 && b[b.length() - i - 1] <= 57) {
                usB += to_string(b[b.length() - i - 1] - '0');
            }
        }

        int ze1 = 0, ze2 = 0;
        for (int i = 0; i < usA.length(); ++i) {
            if (ze1 == i && usA[i] == '0')
                ze1++;
            else
                break;
        }
        for (int i = 0; i < usB.length(); ++i) {
            if (ze2 == i && usB[i] == '0')
                ze2++;
            else
                break;
        }
        string res = karatsubaMul(usA, usB, ze1 + ze2);
        string sigDo = judSignDot(a, b);
        string sig, dt;
        istringstream is(sigDo);
        is >> sig >> dt;

        if (dt != "0")
            res = res.insert(res.length() - atoi(dt.c_str()), 1, '.');

        long firInt = 0;
        long firDot = 0;
        long lastDot = 0;
        long lastInt = 0;
        for (int i = 0; i < res.length(); ++i) {
            if (res[i] != '0') {
                if (res[i] == '.')
                    firDot = i;
                else
                    firInt = i;
                break;
            }
        }
        if (firDot != 0)
            res = res.substr(firDot - 1, res.length());
        else if (firInt != 0)
            res = res.substr(firInt, res.length());

        for (int i = 0; i < res.length(); ++i) {
            if (res[i] == '.') {
                lastDot = i;
                continue;
            }
            if (lastDot != 0 && res[i] != '0')
                lastInt = i;
        }
        if (lastInt > lastDot)
            res = res.substr(0, lastInt + 1);
        else if (lastInt < lastDot)
            res = res.substr(0, lastDot);
        if (sig == "-")
            res = sig + res;
        cout << "Calculation successful! \n";
        cout << "The answer is:" << res << endl;
        return true;
    }
    return false;
}

string karatsubaMul(string n, string m, int zenum) {

    //递归分制的终止条件
    if (n.length() < 4 || m.length() < 4) {
        int numN[n.length()], numM[m.length()];

        int indN = n.length() - 1;
        int indM = m.length() - 1;

        for (int i = n.length() - 1; i >= 0; --i) {
            numN[indN--] = n[n.length() - i - 1] - '0';
        }

        for (int i = m.length() - 1; i >= 0; --i) {
            numM[indM--] = m[m.length() - i - 1] - '0';
        }
        return intMul(numN, numM, n.length(), m.length(), "+ 0");
    }

    // 计算拆分长度
    int sizeN = n.length();
    int sizeM = m.length();
    int half;
    if (sizeN >= sizeM)
        half = (sizeM) / 2;
    else half = (sizeN) / 2;

    //分为a, b, c, d
    string a = n.substr(0, sizeN - half);
    string b = n.substr(sizeN - half, sizeN);
    string c = m.substr(0, sizeM - half);
    string d = m.substr(sizeM - half, sizeM);

    //使用递归
    string p2 = karatsubaMul(a, c, 0);
    string p0 = karatsubaMul(b, d, 0);
    string p1 = myPlus(karatsubaMul(myPlus(a, b, "+"), myPlus(c, d, "+"), 0), myPlus(p0, p2, "+"), "-");

    for (int i = 0; i < half * 2; ++i) {
        p2 += "0";
    }
    for (int i = 0; i < half; ++i) {
        p1 += "0";
    }
    string mypl = myPlus(myPlus(p1, p2, "+"), p0, "+");

    for (int i = 0; i < zenum; i++) {
        mypl = "0" + mypl;
    }
    return mypl;
}

string myPlus(string a, string b, string sign) {
    int maxLen = 0;
    if (a.length() >= b.length())
        maxLen = a.length();
    else
        maxLen = b.length();

    int inA[maxLen];
    int inB[maxLen];
    int ina = a.length() - 1;
    int inb = b.length() - 1;
    for (int i = maxLen - 1; i >= 0; --i) {
        if (ina >= 0) {
            inA[i] = a[ina--] - '0';
        } else
            inA[i] = 0;

        if (inb >= 0) {
            inB[i] = b[inb--] - '0';
        } else
            inB[i] = 0;
    }

    ina = maxLen - 1;
    inb = maxLen - 1;
    int rs[maxLen + 1];
    for (int i = 0; i < maxLen + 1; ++i) {
        rs[i] = 0;
    }
    if (sign == "+") {
        for (int i = maxLen; i >= 1; --i) {
            int temp = rs[i] + inA[ina--] + inB[inb--];
            rs[i] = temp % 10;
            rs[i - 1] = temp / 10;
        }
    } else {
        for (int i = maxLen; i >= 1; --i) {
            if (inA[ina--] + rs[i] >= inB[inb--])
                rs[i] += inA[ina + 1] - inB[inb + 1];
            else {
                rs[i - 1]--;
                rs[i] += inA[ina + 1] + 10 - inB[inb + 1];
            }
        }

    }
    string reu;
    int numb = 0;
    for (int i = 0; i < maxLen + 1; ++i) {
        if (i == numb && rs[i] == 0 && numb != maxLen) {
            numb++;
            continue;
        }

        reu += to_string(rs[i]);
    }
    return reu;
}