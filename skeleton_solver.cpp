#include <iostream>
#include <cstring>
#include <map>
#include <set>
#include <vector>
#include <queue>
#include <utility>
#include <algorithm>
#include <bitset>
#include <cassert>
#include <ctime>

using namespace std;

#define SZN 17
#define x first
#define y second
#define _ ios_base::sync_with_stdio(0); 
// #define debug

typedef pair<int, int> PII;


struct node {
    int mask; // bitmask for state
    int nxt[SZN]; // edges;
} nodes[SZN<<1];

int n, m, orig_n, orig_m;
int mat[SZN][SZN], mat_cpy[SZN][SZN], mat_sol[SZN][SZN];
int M[SZN]; // used to remove duplicates
int orig_row_idx[SZN], orig_col_idx[SZN];
vector<int> duplicate_rows[3][SZN]; // old, new, absorb

// two state perfect phylogeny
int M_mask[SZN<<1]; // sorted col bitmask
int mat2[SZN][SZN<<1];
int num[SZN<<1];
vector<PII> v;
int num_nodes;
bool used[SZN<<1];

// enumerate all children
int addition[SZN];
vector<int> possible_ancestor[SZN];
set<string> solutions;

// prufer
vector<int> edges[SZN];
bool visited[SZN];
int pseq[SZN];
int deg[SZN];
queue<int> q;
int mask[SZN];


bool add_node(int row, int m) {
    int cur = 0;
    for (int i = m-1; i >= 0; --i) {
        if (M_mask[row] & (1<<i)) {
            if (nodes[cur].nxt[i] == -1) {
                if (used[i]) return 0;
                used[i] = 1;
                nodes[cur].nxt[i] = ++num_nodes;
                nodes[num_nodes].mask = nodes[cur].mask | (1<<i);
                memset(nodes[num_nodes].nxt, -1, sizeof(nodes[num_nodes].nxt));
            }
            cur = nodes[cur].nxt[i];
        }
    }

    return 1;
}

void clear() {
    nodes[0].mask = 0;
    memset(nodes[0].nxt, -1, sizeof(nodes[0].nxt));
    memset(used, 0, sizeof(used));
    num_nodes = 0;
}

bool cmp(const PII& a, const PII& b) {
    if (a.x == b.x) return a.y < b.y;
    return a.x > b.x;
}

void add_solution() {
    string sol = "";
    for (int i = 0; i < orig_n; ++i) {
        for (int j = 0; j < orig_m; ++j) {
            sol += (mat_sol[i][j] + '0');
            sol += " ";
        }
        sol += "\n";
    }
    solutions.insert(sol); // are solutions distinct?
}


bool solve_two_state_phylogeny() {
    memcpy(mat_sol, mat_cpy, sizeof(mat_cpy));
    memset(mat2, 0, sizeof(mat2));
    int lose;
    for (int i = 0; i < n; ++i) {
        lose = M[i] ^ mask[possible_ancestor[i][addition[i]]];
        for (int j = 0; j < m; ++j) {
            if (M[i] & (1<<j)) {
                mat2[i][j<<1|1] = 1;
            } else if (lose & (1<<j)) {
                mat2[i][j<<1] = mat2[i][j<<1|1] = 1;
                for (int k : duplicate_rows[0][orig_row_idx[i]])
                    mat_sol[k][orig_col_idx[m-1-j]] = 2;
            }
        }
    }
    memset(num, 0, sizeof(num));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m<<1; ++j)
            num[j] += mat2[i][j];

    for (int j = 0; j < m<<1; ++j)
        v.emplace_back(num[j], j);
    sort(v.begin(), v.end(), cmp);
    for (int i = 0; i < n; ++i) {
        M_mask[i] = 0;
        for (int j = 0; j < v.size(); ++j)
            M_mask[i] = (M_mask[i]<<1) | mat2[i][v[j].y];
    }
    v.clear();

    clear();
    bool good = 1;
    for (int i = 0; i < n; ++i) {
        good = add_node(i, m<<1);
        if (!good) break;
    }

    return good;
}

bool enumerate_add() {
    int cur = 0;
    while (cur < n) {
        addition[cur]++;
        if (addition[cur] >= possible_ancestor[cur].size()) addition[cur++] = 0;
        else break;
    }
    return cur < n;
}

void try_all_combo() {
    memset(addition, 0, sizeof(addition));
    do {
        if (solve_two_state_phylogeny()) 
            add_solution();
    } while (enumerate_add());
}

void test_prufer() {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= m; ++j)
            if ((M[i] & mask[j]) == M[i])
                possible_ancestor[i].push_back(j);
    bool good = 1;
    for (int i = 0; i < n; ++i)
        if (!possible_ancestor[i].size())
            good = 0;
    if (good) try_all_combo();

    for (int i = 0; i < n; ++i)
        possible_ancestor[i].clear();
}

void print_all_completions() {
    for (string s : solutions)
        cout << s << "\n";
    solutions.clear();
}

bool prufer_add() {
    int cur = 0;
    while (cur < m-1) {
        pseq[cur]++;
        if (pseq[cur] > m) pseq[cur++] = 0;
        else break;
    }
    return cur < m-1;
}

void create_prufer() {
    for (int i = 0; i <= m; ++i) deg[i] = 1;
    for (int i = 0; i < m-1; ++i) deg[pseq[i]]++;
    for (int i = 0; i < m-1; ++i) {
        for (int j = 0; j <= m; ++j) {
            if (deg[j] == 1) {
                edges[pseq[i]].push_back(j);
                edges[j].push_back(pseq[i]);
                deg[j] = 0;
                deg[pseq[i]]--;
                break;
            }
        }
    }
    int u = -1;
    for (int i = 0; i <= m; ++i) {
        if (deg[i] == 1) {
            if (u == -1) u = i;
            else {
                edges[u].push_back(i);
                edges[i].push_back(u);
                break;
            }
        }
    }

    memset(visited, 0, sizeof(visited));
    q.push(m);
    visited[m] = 1;
    mask[m] = 0;
    while (q.size()) {
        u = q.front();
        q.pop();
        for (int nxt : edges[u]) {
            if (!visited[nxt]) {
                visited[nxt] = 1;
                mask[nxt] = mask[u] | (1<<nxt);
                q.push(nxt);
            }
        }
    }
}

void delete_prufer() {
    for (int i = 0; i <= m; ++i) 
        edges[i].clear();
}

void solve() {
    memset(pseq, 0, sizeof(pseq));
    do {
        create_prufer();
        test_prufer();
        delete_prufer();
    } while (prufer_add());

    print_all_completions();
}

void remove_cols() {
    memset(num, 0, sizeof(num));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            num[j] += mat[i][j];
    int mm = m;
    for (int j = 0; j < mm; ++j) {
        if (num[j] < 2) {
            --m;
            orig_col_idx[j] = orig_col_idx[m];
            for (int i = 0; i < n; ++i)
                mat[i][j] = mat[i][m];
        }
    }
}

bool remove_rows() {
    map<int, int> rows;
    for (int i = 0; i < n; ++i) {
        M[i] = 0;
        for (int j = 0; j < m; ++j)
            M[i] = M[i]<<1|mat[i][j];
        if (!rows.count(M[i])) rows[M[i]] = orig_row_idx[i];
        else duplicate_rows[2][rows[M[i]]].push_back(orig_row_idx[i]);
    }

    if (rows.size() == n && !rows.count(0)) return 0;
    else {
        rows.erase(0);
        n = 0;
        for (auto it : rows) {
            for (int j = 0; j < m; ++j)
                mat[n][j] = (it.x>>(m-1-j)) & 1;
            duplicate_rows[1][n] = duplicate_rows[0][it.y];
            for (int j : duplicate_rows[2][it.y])
                for (int k : duplicate_rows[0][j])
                    duplicate_rows[1][n].push_back(k);
            orig_row_idx[n++] = it.y;
        }
        for (int i = 0; i < n; ++i) {
            duplicate_rows[0][i] = duplicate_rows[1][i];
            duplicate_rows[1][i].clear();
            duplicate_rows[2][i].clear();
        }
    }

    return 1;
}

void remove_trivial_rows_cols() {
    // want orig_row_idx[i'] = original i before modification
    while (1) {
        remove_cols();
        if (!remove_rows()) break;
    }
}

/* Compile with 
 g++ -std=c++11 -O2 -Wall skeleton_solver.cpp 
 Run with -time flag to check execution time
 Run with -remove_trivial flag to remove trivial rows and columns (don't use this)
 */
int main(int argc, char* argv[]) { _ // disable sync with stdio
    clock_t t = 0;
    if (argc == 2 && !strcmp(argv[1], "-time")) {
        t = clock();
    }
    while (cin >> n >> m) {
        orig_n = n, orig_m = m;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                cin >> mat[i][j];
        memcpy(mat_cpy, mat, sizeof(mat));

        for (int i = 0; i < n; ++i) {
            orig_row_idx[i] = i;
            duplicate_rows[0][i].push_back(i);
        }
        for (int j = 0; j < m; ++j)
            orig_col_idx[j] = j;

        if (argc == 3 && !strcmp(argv[2], "-remove_trivial")) {
            // check for trivial columns and rows and remove them
            remove_trivial_rows_cols();
            if (!n) { // removed all rows
                // only finds one solution if we don't need loss edges
                memcpy(mat_sol, mat_cpy, sizeof(mat_cpy));
                add_solution();
            }
        }

        for (int i = 0; i < n; ++i) {
            M[i] = 0;
            for (int j = 0; j < m; ++j)
                M[i] = M[i]<<1|mat[i][j];
        }

        solve();
    }
    if (argc == 2 && !strcmp(argv[1], "-time")) {
        t = clock() - t;
        cout << "Time elapsed: " << ((double) t / CLOCKS_PER_SEC) << "\n";
    }

    return 0;
}
