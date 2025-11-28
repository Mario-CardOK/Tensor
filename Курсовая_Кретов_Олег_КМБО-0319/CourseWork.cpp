
#include <iostream>
#include <math.h>
#include <vector>
#include <clocale>
#include <fstream>
#include <io.h>
#include <sys/stat.h>
using namespace std;
unsigned long factorial(unsigned long number)
{
    if (number <= 1)
        return 1;
    else
        return number * factorial(number - 1);
}
vector<vector<int>> transpos(vector<int> index)
{
    vector<vector<int>> res;
    vector<int> elem(2);
    int temp;
    int tempJ = -1;
    for (int i = 0; i < index.size() - 1; i++)
    {
        elem[0] = (i + 1);
        for (int j = 0; j < index.size(); j++)
        {
            if (i + 1 == index[j]) { elem[1] = (j + 1); tempJ = j; }
        }
        if (index[i] == index[tempJ]) continue;
        temp = index[i];
        index[i] = index[tempJ];
        index[tempJ] = temp;

        res.push_back(elem);

    }
    return res;
}
int trans_parity(vector<int> trans)
{
    int invers = 0;
    for (int i = 0; i < trans.size(); i++)
    {
        for (int j = i; j < trans.size(); j++)
        {
            if (trans[i] > trans[j] and i < j) invers++;
        }
    }
    if (invers % 2 == 1) return -1;
    else return 1;
}
class Index
{
private:
    vector<int> ind;
    int rang;
    int dim;
    int iter;
    int upperIndex;
public:
    Index(int d, int rg, int upp)
    {
        for (int i = 0; i < rg; i++)
        {
            ind.push_back(1);
        }
        dim = d;
        rang = rg;
        iter = 0;
        upperIndex = upp;
    }
    Index begin()
    {
        Index i(dim,rang, upperIndex);
        return i;
    }
    void operator++(int)
    {
        
        if (iter == pow(dim, rang) - 1)
        {
            cout << "no next item\n";
        }else
            for (int i = (ind.size() - 1); i != -1; i--)
        {
            
            
            if (ind[i] == dim and i != 0)
            {
                ind[i] = 1;
            }
            else
            {
                
                ind[i] = ind[i] + 1;
                iter += 1;
                break;
            }
        }
    }
    vector<int> getElem() { return ind; }
    int IterIn(){ return iter;}
    Index operator+(Index a)
    {
        Index result(dim, rang + a.rang, upperIndex + a.upperIndex);
        vector<int> res;
        for (int i = 0; i < upperIndex; i++)
        {
            res.push_back(ind[i]);
        }
        for (int i = 0; i < a.upperIndex; i++)
        {
            res.push_back(a.ind[i]);
        }
        for (int i = upperIndex; i < rang ; i++)
        {
            res.push_back(ind[i]);
        }
        for (int i = a.upperIndex; i < a.rang; i++)
        {
            res.push_back(a.ind[i]);
        }
        result.ind = res;
        return result;
    }
    int operator[](int i)
    {
        return ind[i];
    }
    Index transp(int a, int b)
    {
        vector<int> TempIndx = ind;
        Index res = *this;
        int a1 = ind[a - 1];
        int b1 = ind[b - 1];
        res.ind[a-1] = b1;
        res.ind[b-1] = a1;
        Index temp(dim, rang, upperIndex);
        for (int i = 0; i < pow(dim, rang); i++)
        {
            if (res.ind == temp.ind)
            {
                res.iter = i;
                break;
            }
            temp++;

        }
        return res;

    }
};

class Tensor
{

private:
    int dimV;
    int rang;
    vector<double long> tensor;
    vector<Index> indexId;
    int upperIndex;
    int lowerIndex;
    int quantityElems;
    Tensor simetrAssist(int posit)// 0 - если нужно симметризация, 1 - антисимметризация
    {
        vector<vector<int>> transLow;//все возможные перестановки нижних индексов
        Tensor res(dimV, 0, rang);
        Tensor temp(dimV, upperIndex, lowerIndex);
        //вектор состоящий из векторов размерности 2
        vector<vector<int>> trLow;
        bool flag;
        if (upperIndex != 0)
        {
            temp = *this;
            while (temp.upperIndex != 0)
            {
                temp = temp.indexLow(1);
            }
        }
        else temp = *this;
        Index indLow(rang, rang, 0);
        Tensor temp2(dimV, 0, rang);
        for (int i = 0; i < pow(rang, rang) - 1; i++)//выбираю только те индексы, в которых не повторов
        {
            flag = 1;
            for (int j = 0; j < rang; j++)
            {
                for (int k = j + 1; k < rang; k++)
                {
                    if (indLow.getElem()[j] == indLow.getElem()[k]) { flag = 0; break; }
                }
            }
            if (flag) transLow.push_back(indLow.getElem());

            indLow++;
        }
        for (int k = 0; k < transLow.size(); k++)
        {
            temp2 = temp;
            trLow = transpos(transLow[k]);
            for (int j = 0; j < trLow.size(); j++)
            {

                temp2 = temp2.transposLow(trLow[j][0], trLow[j][1]);
            }
            if (posit == 0) res += temp2;
            else res += temp2 * trans_parity(transLow[k]);

        }
        double long  e = factorial(rang);
        res = res * (1 / e);
        return res;
    }

public:
    Tensor(int dim, int uIndex, int lIndex)
    {
        if (dim == 0) { cerr << "Размерность не может быть равной 0\n"; exit(-1); }
        if(uIndex < 0 or lIndex < 0 or dim < 1) { cerr << "Проблемы ввода характеристик тензора\n"; exit(-1); }
        dimV = (dim);
        upperIndex = (uIndex);
        lowerIndex = (lIndex);
        rang = upperIndex + lowerIndex;
        quantityElems = pow(dimV, rang);
        tensor = vector<double long> (quantityElems, 0);

        Index index(dimV, rang, upperIndex);
        indexId.push_back(index);
        for (int i = 0; i < quantityElems-1; i++)
        {
            index++;
            indexId.push_back(index);  
        }
    }
    void addElems(vector<long double> coord)
    {
        Index index(dimV, rang, upperIndex);
        for (int i = 0; i < coord.size(); i++) tensor[i] = coord[i];
    }
    double long getElem(vector<int> comp)
    {
        
        for (int i = 0; i < quantityElems; i++)
        {
            if (indexId[i].getElem() == comp) return tensor[i];
        }
    }
    void setElem(vector<int>  comp, double item)
    {
        for (int i = 0; i < quantityElems; i++)
        {
            if (comp == indexId[i].getElem()) { tensor[i] = item; break; }
        }
    }
    void setElem(Index comp, double item)
    {
        for (int i = 0; i < quantityElems; i++)
        {
            if (comp.getElem() == indexId[i].getElem()) { tensor[i] = item; break; }
        }
    }
    double long operator[](int i)
    {
        return tensor[i];
    }
    Tensor operator+(Tensor t)
    {
        for (int i = 0; i < quantityElems; i++)
        {
            tensor[i] += t[i];
        }
        return *this;
    }
    Tensor operator-(Tensor t)
    {
        for (int i = 0; i < quantityElems; i++)
        {
            tensor[i] -= t[i];
        }
        return *this;
    }
    Tensor operator+=(Tensor t)
    {
        for (int i = 0; i < quantityElems; i++)
        {
            tensor[i] += t[i];
        }
        return *this;
    }
    Tensor operator*=(double long l)
    {
        for (int i = 0; i < quantityElems; i++)
        {
            tensor[i] *= l;
        }
        return *this;
    }
    Tensor operator*(double long l)
    {
        Tensor res(dimV, upperIndex, lowerIndex);
        for (int i = 0; i < quantityElems; i++)
        {
            res.tensor[i] = tensor[i] * l;
        }
        return res;
    }
    Tensor operator*(Tensor tens)
    {
        int newRang = rang + tens.rang;
        int newUpIndex = upperIndex + tens.upperIndex;
        int lowIndex = lowerIndex + tens.lowerIndex;
        int quantity = quantityElems * tens.quantityElems;
        vector<Index> newIndexId;
        Index index(dimV, newRang, newUpIndex);
        newIndexId.push_back(index);
        for (int i = 0; i < quantity-1; i++)
        {
            index++;
            newIndexId.push_back(index);
        }
        Tensor res(dimV, newUpIndex, lowIndex);
        res.indexId = newIndexId;
        res.quantityElems = quantity;
        res.rang = newRang;
        vector<double long> v(quantity, 0);
        res.tensor = v;
        double elem;
        Index comp(dimV,newRang, newUpIndex);

        for (int i = 0; i < quantityElems; i++)
        {
            for (int j = 0; j < tens.quantityElems; j++)
            {
                elem = tensor[i] * tens.tensor[j];
                comp = indexId[i] + tens.indexId[j];
                res.setElem(comp, elem);
            }
        }

        return res;
    }
    Tensor operator*=(Tensor tens)
    {
        int newRang = rang + tens.rang;
        int upIndex = upperIndex + tens.upperIndex;
        int lowIndex = lowerIndex + tens.lowerIndex;
        int quantity = quantityElems * tens.quantityElems;
        vector<Index> newIndexId;
        Index index(dimV, newRang, upIndex);

        for (int i = 0; i < quantity; i++)
        {
            newIndexId.push_back(index);
            index++;
        }
        Tensor res(dimV, upIndex, lowIndex);
        res.indexId = newIndexId;
        res.quantityElems = quantity;
        res.rang = newRang;
        vector<double long> v(quantity, 0);
        res.tensor = v;


        double elem;
        vector<int> comp(newRang, 6);
        int k = 0;
        for (int i = 0; i < quantityElems; i++)
        {
            for (int j = 0; j < quantityElems; j++)
            {
                elem = tensor[i] * tens.tensor[j];
                for (int l = 0; l < newRang; l++)
                {
                    if (l < upperIndex) { comp[l] = indexId[i].getElem()[l % rang]; }
                    else if ((l < (upperIndex + tens.upperIndex))) { comp[l] = tens.indexId[j].getElem()[l % tens.rang]; }
                    else if ((l < (upperIndex + tens.upperIndex + lowerIndex))) { comp[l] = indexId[i].getElem()[l % rang]; }
                    else if (l < upperIndex + tens.upperIndex + lowerIndex + tens.lowerIndex) { comp[l] = tens.indexId[j].getElem()[l % tens.rang]; }

                }

                res.setElem(comp, elem);
                k++;
            }
        }
        *this = res;
        return *this;
    }
    Tensor contraction(int m, int n)
    {
        if (m > upperIndex or n > lowerIndex)
        {
            cerr << "Index error";
        }
        int newRang = rang - 2;
        int newUpper = upperIndex - 1;
        int newLower = lowerIndex - 1;

        Tensor res(dimV, newUpper, newLower);
        Tensor elem(dimV, newUpper, newLower);
        Index newIndx(dimV, newRang, newUpper);

        for (int i = 0; i < dimV; i++)
        {
            newIndx = newIndx.begin();
            for (int j = 0; j < quantityElems; j++)
            {
                if (indexId[j].getElem()[m-1] == i+1 and indexId[j].getElem()[upperIndex +n - 1] == i+1)
                {
                    elem.setElem(newIndx.getElem(), tensor[indexId[j].IterIn()]);
                    if(newIndx.IterIn() != res.quantityElems - 1) newIndx++;
                }
            }
            res += elem;
        }
        return res;
    }
    Tensor indexUp(int ind)
    {
        Tensor res(dimV, upperIndex+1, lowerIndex-1);
        Tensor metric(dimV, 2, 0);
        vector<double long> MetTens(dimV*dimV);
        for (int i = 0; i < MetTens.size(); i+= (dimV + 1)) MetTens[i] = 1;
        
        metric.addElems(MetTens);
        res = metric * (*this);
        return res.contraction(1,ind);
    }
    Tensor indexLow(int ind)
    {
        Tensor res(dimV, upperIndex - 1, lowerIndex + 1);
        Tensor metric(dimV, 0, 2);
        vector<double long> MetTens(dimV * dimV);
        for (int i = 0; i < MetTens.size(); i += (dimV + 1)) MetTens[i] = 1;
        metric.addElems(MetTens);
        res = metric * (*this);
        return res.contraction(ind, 1);
    }
    Tensor transposUp(int a, int b)
    {
        Tensor res = *this;
        Index tempInd(dimV, rang, upperIndex);
        for (int i = 0; i < quantityElems; i++)
        {
            tempInd = indexId[i].transp(a,b);
            res.tensor[i] = tensor[tempInd.IterIn()];
        }
        return res;
    }
    Tensor transposLow(int a, int b)
    {
        Tensor temp = *this;
        Index tempInd(dimV, rang, upperIndex);
        for (int i = 0; i < quantityElems; i++)
        {
            tempInd = indexId[i].transp(upperIndex + a, upperIndex + b);
            temp.tensor[i] = tensor[tempInd.IterIn()];
        }
        return temp;
    }
    Tensor simetrisation()
    { 
        return (*this).simetrAssist(0);
    }
    Tensor antisimetr()
    {
        
        return (*this).simetrAssist(1);
    }
    int Rang() { return rang; }
    int UperIndx() { return upperIndex; }
    int LowIndx() { return lowerIndex; }
    friend ostream& operator<< (std::ostream& out, const Tensor& tensor);
    friend istream& operator>> (std::istream& in, Tensor& point);
    friend ifstream& operator>> (ifstream& f, Tensor& T);
    friend ofstream& operator<<(ofstream& out, const Tensor& T);

};

ostream& operator<< (ostream& out, const Tensor& T)
{
    for (int i = 0; i < T.tensor.size(); i++)
    {
        out << T.tensor[i] << "  ";
        if (T.upperIndex == 1 and T.lowerIndex == 0) out << endl;
        else if ((i + 1)% T.dimV == 0) out << endl;
    }
    return out;
}
istream& operator>> (istream& in, Tensor& T)
{
    cout << "Введите размерность пространства, кол-во вверхних и нижних индексов\n";
    try
    {
        int dim, up, low = 0;
        in >> dim;
        in >> up;
        in >> low;
        T = Tensor(dim, up, low);
        cout << "Теперь вводите значения:\n";
        for (int i = 0; i < T.tensor.size(); i++)
        {
            in >> T.tensor[i];
        }
    }
    catch (...)
    {
        cerr << "Ошибка ввода, перезапустите программу и попробуйте еще раз"; exit(-1);
    }
    
    return in;
}
ifstream& operator>> (ifstream& f, Tensor& T)
{
    try
    {
        int dim, up, low = 0;
        f >> dim;
        f >> up;
        f >> low;
        T = Tensor(dim, up, low);
        for (int i = 0; i < T.tensor.size(); i++)
        {
            f >> T.tensor[i];
        }
    }
    catch (...)
    {
        cerr << "Ошибка ввода, перезапустите программу и попробуйте еще раз\n"; exit(-1);
    }
    
    return f;
}
ofstream& operator<<(ofstream& out, const Tensor& T)
{
    for (int i = 0; i < T.tensor.size(); i++)
    {
        out << T.tensor[i] << "  ";
        if(T.upperIndex == 1 and T.lowerIndex == 0 ) out << endl;
        else if ((i + 1) % T.dimV == 0) out << endl;
    }
    return out;
}

int main()
{
    
    setlocale(LC_ALL, "Russian");
    cout << "*Исходный тензор при выполнении операций остается неизменным.\n";
    cout << "Выберите, как вы хотите получить тензор:\nВвести самому - 0    Загрузить из input.txt - 1    Демонстрация - 2\n";
    ofstream of("output.txt", ios::trunc);
    int flag;
    cin >> flag;
    Tensor t(1,0,0);
    
    while (true)
    {
        if (flag == 0)
        {
            cin >> t;
            cout << "Вы ввели:\n";
            cout << t << endl;
            break;
            
        }
        else if (flag == 1)
        {
            
            
            ifstream f("input.txt");
            struct stat sb;
            if (stat("input.txt", &sb) == -1 && errno == ENOENT) { cerr << "Не удалось открыть текстовый документ input.txt\n"; exit(-1); }
            if (f.peek() == EOF) { cerr << "Пустой input.txt\n"; exit(-1); }
            f >> t;
            cout << "Было считано:\n";
            cout << t << endl;

            break;
            
        }
        else if (flag == 2)
        {
            vector<double long> v = { 1, 2 , 3, 4 };
            t = Tensor(2, 1, 1);
            t.addElems(v);
            cout << "Демонстрационный тензор:\nРазмерность - 2 , тип (1,1)\n";
            cout << "Его значения:\n";
            cout <<  t << endl;
            break;
        }
        else
        {
            cout << "Что-то пошло не так, попробуйте еще раз\n";
            cout << "Выберите, как вы хотите получить тензор:\nВвести самому - 0    Загрузить из input.txt - 1    Демонстрация - 2\n";
            cin >> flag;
        }
    }
    
    
    
    cout << "Тензор умноженый на 10:\n";
    cout << t*10 << endl;
    of << "Тензор умноженый на 10:\n";
    of << t * 10 << endl;

    if (t.Rang() > 1)
    {
        cout << "Симметричный тензор:\n";
        cout << t.simetrisation() << endl;
        of << "Симметричный тензор:\n";
        of << endl << t.simetrisation();
        cout << "Антисимметричный тензор:\n";
        cout << t.antisimetr() << endl;
        of << "Антисимметричный тензор:\n";
        of << endl << t.antisimetr();
        cout << "Симметричный тензор + Антисимметричный тензор:\n";
        cout << t.antisimetr() + t.simetrisation() << endl;
        of << "Симметричный тензор + Антисимметричный тензор:\n";
        of << endl << t.antisimetr() + t.simetrisation();
    }
    if (t.UperIndx() > 0 and t.LowIndx() > 0)
    {
        cout << "Свертка:\n";
        cout << t.contraction(1, 1) << endl;
        of << "Свертка:\n";
        of << endl << t.contraction(1, 1);
    }
    if (t.UperIndx() > 0)
    {
        cout << "Опускание первого верхнего индекса:\n";
        cout << t.indexLow(1);
        of << "Опускание первого верхнего индекса:\n";
        of << endl << t.indexLow(1);
    }
    if (t.LowIndx() > 0)
    {
        cout << "Поднятие первого нижнего индекса:\n";
        cout << t.indexUp(1);
        of << "Поднятие первого нижнего индекса:\n";
        of << endl << t.indexUp(1);
    }
    if (t.UperIndx() > 1)
    {
        cout << "Перестановка 1 и 2 верхних индексов:\n";
        cout << t.transposUp(1, 2);
        of << "Перестановка 1 и 2 верхних индексов:\n";
        of << endl << t.transposUp(1, 2);
    }
    if (t.LowIndx() > 1)
    {
        cout << "Перестановка 1 и 2 нижних индексов:\n";
        cout << t.transposLow(1, 2);
        of << "Перестановка 1 и 2 нижних индексов:\n";
        of << endl << t.transposLow(1, 2);
    }

    cout << "Тензор умноженный сам на себя:\n";
    cout << t * t << endl;
    of << "\nТензор умноженный сам на себя:\n";
    of << t * t << endl;
    char c;
    cin >> c;
    return 0;
}