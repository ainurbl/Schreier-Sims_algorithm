#include <bits/stdc++.h>

using namespace std;

int lcm(int ret, int x) { return ret / __gcd(ret, x) * x; }

class Permutation {
 public:
  int n{};
  vector<int> perm;  // perm[0] is fictive, permutation represented as table
                     // 1-indexed
  explicit Permutation(int n) {
    this->n = n;
    perm.resize(n + 1);
    for (int i = 1; i <= n; i++) {
      perm[i] = i;
    }
  }

  Permutation() = default;

  int size() const { return (int)perm.size() - 1; }

  int &operator[](int id) { return perm[id]; }

  Permutation operator*(
      Permutation other) {  // multiplication of two permutations
    Permutation result(n);
    for (int i = 1; i <= n; i++) {
      result[i] = perm[other[i]];
    }
    return result;
  }

  int operator*(int base) { return perm[base]; }

  void
  print() {  // printing permutation in "good" format (multiplication of cycles)
    vector<int> used(size() + 1, 0);
    bool isId = true;
    for (int i = 1; i <= size(); i++) {
      if (used[i]) continue;
      int cur = i;
      vector<int> cycle;
      cycle.push_back(cur);
      used[cur] = 1;
      cur = perm[cur];
      while (cur != i) {
        cycle.push_back(cur);
        used[cur] = 1;
        cur = perm[cur];
      }
      if (cycle.size() > 1) {  // cycle is not empty
        isId = false;
        cout << "(";
        for (int j = 0; j < cycle.size() - 1; j++) {
          cout << cycle[j] << ", ";
        }
        cout << cycle.back() << ") ";
      }
    }
    if (isId) {
      cout << "id ";
    }
    cout << "\n";
  }

  vector<int> cyclicType() {
    vector<int> res;  // sizes of cicles
    vector<int> used(size() + 1, 0);
    bool isId = true;
    for (int i = 1; i <= size(); i++) {
      if (used[i]) continue;
      int cur = i;
      vector<int> block;
      block.push_back(cur);
      used[cur] = 1;
      cur = perm[cur];
      while (cur != i) {
        block.push_back(cur);
        used[cur] = 1;
        cur = perm[cur];
      }
      if (block.size() > 1) {
        res.push_back(block.size());
      }
    }
    return res;
  }

  bool isNeutral() {
    for (int i = 1; i <= n; i++) {
      if (perm[i] != i) return false;
    }
    return true;
  }

  int getOrd() {
    vector<int> type = cyclicType();
    if (type.empty()) {
      return 1;
    }
    int ret = type[0];
    for (auto &x : type) {
      ret = lcm(ret, x);
    }
    return ret;
  }

  Permutation inversed() {  // permutation inversion
    return binPower(getOrd() - 1);
  }

  Permutation binPower(int power) {
    Permutation res(size());
    Permutation val = *this;
    while (power) {
      if (power % 2 == 1) {
        res = res * val;
      }
      power /= 2;
      val = val * val;
    }
    return res;
  }

  bool operator<(const Permutation &other) const { return perm < other.perm; }
};

class SchreiersTree {
 public:
  int root{}, n{};
  vector<pair<int, int>> tree;  // tree[i] - set of pairs <vertex v, index j
                                // such that v = perm[j] * i>

  vector<Permutation> perm;
  vector<Permutation> pathFromRoot;

  const pair<int, int> EMPTY = {-1, -1};

  SchreiersTree() = default;

  SchreiersTree(vector<Permutation> generatingSet, int base,
                int n) {  // building Schreier's tree
    root = base;
    this->n = n;
    tree.resize(n + 1, EMPTY);
    pathFromRoot.resize(n + 1, Permutation(n));
    perm = move(generatingSet);
    build();
  }

  void build() {
    queue<int> order;
    order.push(root);
    vector<bool> used(n + 1);
    used[root] = true;
    while (!order.empty()) {
      int v = order.front();
      order.pop();
      for (int i = 0; i < perm.size(); i++) {
        int newV = perm[i] * v;
        if (used[newV]) continue;
        pathFromRoot[newV] = perm[i] * pathFromRoot[v];
        assert(pathFromRoot[newV] * root == newV);
        order.push(newV);
        used[newV] = true;
        assert(tree[newV] == EMPTY);
        tree[newV] = {v, i};
      }
    }
  }

  Permutation getPerm(int val) {
    if (tree[val] == EMPTY) return Permutation(n);
    return perm[tree[val].second];
  }

  void test() {
    for (int i = 1; i <= n; i++) {
      if (tree[i] == EMPTY) continue;
      int val = i;
      Permutation sigma = perm[tree[val].second];
      val = tree[val].first;
      while (val != root) {
        sigma = sigma * perm[tree[val].second];
        val = tree[val].first;
      }
      assert(sigma.inversed() * i == root);
    }
  }

  vector<int> getOrbit() {
    vector<int> res;
    res.push_back(root);
    for (int i = 1; i <= n; i++) {
      if (tree[i] != EMPTY) res.push_back(i);
    }
    return res;
  }

  bool isElementInOrbit(int a) {
    vector<int> orbit = getOrbit();
    return find(orbit.begin(), orbit.end(), a) != orbit.end();
  }

  set<Permutation> schreiersLemma() {
    set<Permutation> generator;
    assert(tree[root] == EMPTY);
    for (int y = 1; y <= n; y++) {
      if (tree[y] == EMPTY && y != root) continue;
      for (auto &s : perm) {
        int sy = s * y;
        Permutation h_sy = pathFromRoot[sy];
        Permutation h_y = pathFromRoot[y];
        Permutation toAdd = h_sy.inversed() * s * h_y;
        if (toAdd.isNeutral()) continue;
        generator.insert(toAdd);
      }
    }
    return generator;
  }
};

template <class T>
vector<T> toVector(set<T> a) {
  vector<T> ret;
  ret.reserve(a.size());
  for (auto &x : a) ret.push_back(x);
  return ret;
}

class FullChain {
 public:
  int n;
  vector<int> base;  // base of set
  vector<vector<Permutation>>
      strongGeneratingSets;  // strong set of generators G are in [0],
                             // [i] - generators of Stab^{G}_{(b_1, ..., b_i)}

  vector<SchreiersTree> trees;  // Schreier's trees

  FullChain(const vector<int> &base,
            vector<Permutation> generatingSet) {  // Schreier-Sims algorithm
    n = generatingSet[0].size();
    this->base = base;
    strongGeneratingSets.push_back(generatingSet);
    for (auto &b_i : base) {
      trees.emplace_back(strongGeneratingSets.back(), b_i, n);
      strongGeneratingSets.push_back(toVector(trees.back().schreiersLemma()));
      normalize(strongGeneratingSets.back());
      assert(strongGeneratingSets.back().size() <= n * n);
      for (auto &perm : strongGeneratingSets.back()) {
        assert(perm[b_i] == b_i);
      }
    }
  }

  FullChain getFullChainFor(int i) {  // getting full chain of G_i
    FullChain ret = *this;
    ret.trees.clear();
    ret.strongGeneratingSets.clear();
    ret.base.clear();
    for (int id = i; id < n; id++) {
      if (trees.size() > id) ret.trees.push_back(trees[id]);
      if (strongGeneratingSets.size() > id)
        ret.strongGeneratingSets.push_back(strongGeneratingSets[id]);
      if (base.size() > id) ret.base.push_back(base[id]);
    }
    return ret;
  }

  vector<int> getOrbitForb1() {  // O(b_1) as set
    return trees[0].getOrbit();
  }

  static bool checkIsItFullChain() { return true; }

  size_t getOrder() {  // order of group G
    size_t res = 1;
    for (auto &x : trees) {
      res *= x.getOrbit().size();
    }
    return res;
  }

  void normalize(vector<Permutation> &generatingSet) {
    set<Permutation> newGeneratingSet;
    vector<vector<Permutation>> used(
        n + 1, vector<Permutation>(n + 1, Permutation(1)));
    for (auto &s : generatingSet) {
      Permutation newS = s;
      for (int i = 0; i < base.size(); i++) {
        int x = newS * base[i];
        if (x == base[i]) continue;
        if (used[base[i]][x].n != 1) {
          newS = newS.inversed() * used[base[i]][x];
        } else {
          used[base[i]][x] = newS;
          newGeneratingSet.insert(newS);
          break;
        }
      }
    }
    generatingSet = toVector(newGeneratingSet);
  }

  bool contains(Permutation a) {  // is a in G
    Permutation val = std::move(a);
    for (int i = 0; i < base.size(); i++) {
      int element = val[base[i]];
      if (!trees[i].isElementInOrbit(element)) {
        return false;
      }
      Permutation b = trees[i].getPerm(element).inversed();
      val = b * val;
    }
    return val.isNeutral();
  }

  vector<Permutation> get(Permutation a) {  // present a using generators of G
    vector<Permutation> ret;
    Permutation val = std::move(a);
    for (int i = 0; i < base.size(); i++) {
      int element = val[base[i]];
      if (!trees[i].isElementInOrbit(element)) {
        return {Permutation(n)};
      }
      Permutation b = trees[i].getPerm(element).inversed();
      val = b * val;
      if (!val.isNeutral()) ret.push_back(val);
    }
    if (ret.empty()) ret.emplace_back(n);
    return ret;
  }
};

class Tests {
 public:
  static void runAll() {
    testTree1();
    testTree2();
    testTree3();
    testChainS4();
    testChainV4();
    testChainA5();
    testChainS5();
    testChainD6();
    testChainCornerCase();
    testIdInS4();
    test123InS4();
    test123InA5();
    test1234NotInA5();
    test123NotInV4();
    test12_34InV4();
  }

 private:
  static void testTree1() {
    Permutation a(8), b(8), c(8), d(8);
    a[1] = 2;
    a[2] = 1;
    a[5] = 6;
    a[6] = 7;
    a[7] = 8;
    a[8] = 5;
    b[2] = 3;
    b[3] = 2;
    c[3] = 4;
    c[4] = 3;
    d[4] = 8;
    d[8] = 4;
    SchreiersTree tree({b, d, a, c}, 1, 8);
    tree.test();
    vector<int> expect = {1, 2, 3, 4, 5, 6, 7, 8};
    assert(tree.getOrbit() == expect);
  }

  static void testTree2() {
    Permutation a(8), b(8), c(8), d(8);
    a[1] = 2;
    a[2] = 1;
    a[5] = 6;
    a[6] = 7;
    a[7] = 8;
    a[8] = 5;
    b[2] = 3;
    b[3] = 2;
    c[3] = 4;
    c[4] = 3;
    d[4] = 3;
    d[3] = 4;
    SchreiersTree tree({b, d, a, c}, 1, 8);
    tree.test();
    vector<int> expect = {1, 2, 3, 4};
    assert(tree.getOrbit() == expect);
  }

  static void testTree3() {
    Permutation a(7), b(7), c(7), d(7);
    a[1] = 2;
    a[2] = 3;
    a[3] = 1;
    b[1] = 4;
    b[4] = 5;
    b[5] = 1;
    c[5] = 7;
    c[7] = 5;
    d[3] = 6;
    d[6] = 3;
    SchreiersTree tree({a, d, c, b}, 1, 7);
    tree.test();
    vector<int> expect = {1, 2, 3, 4, 5, 6, 7};
    assert(tree.getOrbit() == expect);
  }

  static void testChainS4() {
    vector<int> base = {1, 2, 3};
    Permutation a(4), b(4), c(4);
    a[2] = 1;
    a[1] = 2;
    b[2] = 3;
    b[3] = 2;
    c[3] = 4;
    c[4] = 3;
    FullChain fullChain(base, {a, b, c});
    cout << "Order of S_4 is " << fullChain.getOrder() << "\n";
    assert(fullChain.getOrder() == 24);
  }

  static void testChainV4() {
    vector<int> base = {1, 2, 3, 4, 5, 6};
    Permutation a(6), b(6);
    a[1] = 2;
    a[2] = 1;
    a[3] = 4;
    a[4] = 3;
    b[1] = 3;
    b[3] = 1;
    b[2] = 4;
    b[4] = 2;
    FullChain fullChain(base, {a, b});
    cout << "Order of V_4 is " << fullChain.getOrder() << "\n";
    assert(fullChain.getOrder() == 4);
  }

  static void testChainA5() {
    vector<int> base = {5, 4, 2, 3, 1};
    Permutation a(5), b(5);
    a[1] = 2;
    a[2] = 3;
    a[3] = 1;
    b[1] = 2;
    b[2] = 3;
    b[3] = 4;
    b[4] = 5;
    b[5] = 1;
    FullChain fullChain(base, {a, b});
    cout << "Order of A_5 is " << fullChain.getOrder() << "\n";
    assert(fullChain.getOrder() == 60);
  }

  static void testChainS5() {
    vector<int> base = {5, 4, 2, 3, 1};
    Permutation a(5), b(5), c(5);
    a[1] = 2;
    a[2] = 3;
    a[3] = 1;
    b[1] = 2;
    b[2] = 3;
    b[3] = 4;
    b[4] = 5;
    b[5] = 1;
    c[1] = 2;
    c[2] = 1;
    FullChain fullChain(base, {c, a, b});
    cout << "Order of S_5 is " << fullChain.getOrder() << "\n";
    assert(fullChain.getOrder() == 120);
  }

  static void testChainD6() {
    vector<int> base = {5, 4, 2, 3, 1, 6};
    Permutation a(6), b(6);
    a[2] = 6;
    a[6] = 2;
    a[3] = 5;
    a[5] = 3;
    b[1] = 2;
    b[2] = 3;
    b[3] = 4;
    b[4] = 5;
    b[5] = 6;
    b[6] = 1;
    FullChain fullChain(base, {b, a});
    cout << "Order of D_6 is " << fullChain.getOrder() << "\n";
    assert(fullChain.getOrder() == 12);
  }

  static void testChainCornerCase() {
    vector<int> base = {5, 4, 2, 3, 1, 6};
    Permutation a(6), b(6);
    FullChain fullChain(base, {b, a});
    cout << "Order of ID is " << fullChain.getOrder() << "\n";
    assert(fullChain.getOrder() == 1);
  }

  static void testIdInS4() {
    vector<int> base = {1, 2, 3};
    Permutation a(4), b(4), c(4);
    a[2] = 1;
    a[1] = 2;
    b[2] = 3;
    b[3] = 2;
    c[3] = 4;
    c[4] = 3;
    FullChain fullChain(base, {a, b, c});
    Permutation id(4);
    cout << "ID is in S_4\n";
    assert(fullChain.contains(id));
  }

  static void test123InS4() {
    vector<int> base = {1, 2, 3};
    Permutation a(4), b(4), c(4);
    a[2] = 1;
    a[1] = 2;
    b[2] = 3;
    b[3] = 2;
    c[3] = 4;
    c[4] = 3;
    FullChain fullChain(base, {a, b, c});
    Permutation sigma(4);
    sigma[1] = 2;
    sigma[2] = 3;
    sigma[3] = 1;
    cout << "(123) is in S_4\n";
    assert(fullChain.contains(sigma));
  }

  static void test123InA5() {
    vector<int> base = {5, 4, 2, 3, 1};
    Permutation a(5), b(5);
    a[1] = 2;
    a[2] = 3;
    a[3] = 1;
    b[1] = 2;
    b[2] = 3;
    b[3] = 4;
    b[4] = 5;
    b[5] = 1;
    FullChain fullChain(base, {a, b});
    Permutation t(5);
    t[1] = 2;
    t[2] = 3;
    t[3] = 1;
    cout << "(123) is in A_5\n";
    assert(fullChain.contains(t));
  }

  static void test1234NotInA5() {
    vector<int> base = {5, 4, 2, 3, 1};
    Permutation a(5), b(5);
    a[1] = 2;
    a[2] = 3;
    a[3] = 1;
    b[1] = 2;
    b[2] = 3;
    b[3] = 4;
    b[4] = 5;
    b[5] = 1;
    FullChain fullChain(base, {a, b});
    Permutation t(5);
    t[1] = 2;
    t[2] = 3;
    t[3] = 4;
    t[4] = 1;
    cout << "(1234) is not in A_5\n";
    assert(!fullChain.contains(t));
  }

  static void test123NotInV4() {
    vector<int> base = {1, 2, 3, 4, 5, 6};
    Permutation a(6), b(6);
    a[1] = 2;
    a[2] = 1;
    a[3] = 4;
    a[4] = 3;
    b[1] = 3;
    b[3] = 1;
    b[2] = 4;
    b[4] = 2;
    FullChain fullChain(base, {a, b});
    Permutation t(6);
    t[1] = 2;
    t[2] = 3;
    t[3] = 1;
    cout << "(123) is not in V_4\n";
    assert(!fullChain.contains(t));
  }

  static void test12_34InV4() {
    vector<int> base = {1, 2, 3, 4, 5, 6};
    Permutation a(6), b(6);
    a[1] = 2;
    a[2] = 1;
    a[3] = 4;
    a[4] = 3;
    b[1] = 3;
    b[3] = 1;
    b[2] = 4;
    b[4] = 2;
    FullChain fullChain(base, {a, b});
    Permutation t(6);
    t[1] = 2;
    t[2] = 1;
    t[3] = 4;
    t[4] = 3;
    cout << "(12)(34) is in V_4\n";
    assert(fullChain.contains(t));
  }
};

class SolveTask {  // task
                   /*
                    * наши перестановки такие: [1, ..., 54] -> [3, ..., 54, 2, 1]
                    *                          [1, ..., 18, 19, ..., 36, 37, ..., 54] ->
                    * перестановка блоков по 18 строим полную цепочку на базе {54, 53...} мы
                    * знаем доп информацию про искомую группу - последний элемент остается на
                    * месте поэтому на самом деле нас интересует такой вопрос: найти все
                    * возможные карты, которые могут быть на предпоследнем месте в Stab_{54}^{G}
                    * где G - подгруппа S_{54}, порожденная нашими перестановками
                    * заметим, что последнее - в точности вопрос о поиске орбиты
                    * так и действуем
                    */
 public:
  Permutation gen(int a, int b, int c) {
    Permutation perm(54);
    for (int i = 1; i <= 18; i++) {
      perm[i] = a + i - 1;
      perm[i + 18] = b + i - 1;
      perm[i + 36] = c + i - 1;
    }
    return perm;
  }

  vector<Permutation> gen() {
    vector<Permutation> ret;
    ret.push_back(gen(1, 37, 19));
    ret.push_back(gen(19, 1, 37));
    ret.push_back(gen(19, 37, 1));
    ret.push_back(gen(37, 1, 19));
    ret.push_back(gen(37, 19, 1));
    return ret;
  }

  SolveTask() {
    Permutation millers(54);
    millers[1] = 54;
    for (int i = 54; i >= 4; i -= 2) {
      millers[i] = i - 2;
    }
    millers[2] = 53;
    for (int i = 53; i >= 3; i -= 2) {
      millers[i] = i - 2;
    }
    vector<Permutation> stierlitz = gen();
    generatingSet = stierlitz;
    generatingSet.push_back(millers);

    cout << "==============Permutations==============\n";

    for (auto &x : generatingSet) {
      x.print();
    }
    cout << "=========================================\n";
    for (int i = 54; i >= 1; i--) {
      base.push_back(i);
    }
  }

  vector<Permutation> generatingSet;
  vector<int> base;

  void showResults() const {
    chrono::time_point<chrono::system_clock> start =
        chrono::system_clock::now();

    FullChain fullChain(base, generatingSet);

    vector<Permutation> Stab_54_generator = fullChain.strongGeneratingSets[1];

    vector<int> orbit_of_53_in_Stab_54 =
        SchreiersTree(Stab_54_generator, 53, 54).getOrbit();
    sort(orbit_of_53_in_Stab_54.begin(), orbit_of_53_in_Stab_54.end());
    cout << "Порядок группы: " << fullChain.getOrder() << "\n";
    for (auto &x : orbit_of_53_in_Stab_54) {
      cout << x << " ";
    }
    cout << "\n";
    cout << "Орбита предпоследней карты в стабилизаторе последней "
            "карты состоит только из предпоследней карты, поэтому "
            "ответ на задачу - Штирлиц однозначно знает карту.\n";
    cout << "Time spent: "
         << chrono::duration_cast<chrono::milliseconds>(
                chrono::system_clock::now() - start)
                .count()
         << " ms\n";
  }
};

class StrongTest {
 public:
  StrongTest() {
    for (int i = 1; i <= 26; i++) base.push_back(i);
    vector<Permutation> a_54(2, Permutation(26));
    a_54[0][1] = 2;
    a_54[0][2] = 3;
    a_54[0][3] = 1;
    for (int i = 2; i <= 25; i++) a_54[1][i] = i + 1;
    a_54[1][26] = 2;
    generatingSet = a_54;
  }

  vector<Permutation> generatingSet;
  vector<int> base;

  void showResults() const {
    chrono::time_point<chrono::system_clock> start =
        chrono::system_clock::now();

    FullChain fullChain(base, generatingSet);

    cout << "A_26 is calculated.\n";

    cout << boolalpha;

    Permutation check(26);
    check[1] = 2;
    check[2] = 1;
    cout << "A_26 contains (12)? " << fullChain.contains(check) << "\n";
    check[12] = 7;
    check[7] = 12;
    cout << "A_26 contains (1 2)(12 7)? " << fullChain.contains(check) << "\n";

    cout << "Time spent: "
         << chrono::duration_cast<chrono::milliseconds>(
                chrono::system_clock::now() - start)
                .count()
         << " ms\n";
  }
};

int main() {
  cout << "==============Tests==============\n";
  Tests::runAll();  // optional
  cout << "=================================\n";
  SolveTask().showResults();
  StrongTest().showResults();  // optional
  return 0;
}
