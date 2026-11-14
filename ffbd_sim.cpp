// ffbd_sim.cpp -- exact simulation for simulation-based calibration (validation step 2).
// Forward birth-death-fossil-sampling from the root (2 lineages at the root age), rejection on the observed configuration:
// n extant sampled tips, m fossils, labelled tips consistent with the fixed extant topology ((A,B),(C,D)) and the clade constraints,
// fossil ages inside their windows. Rejection is joint with the rate draws, so accepted draws are exact samples from
//   p(lambda, mu, psi, root, rates, m, tree | n, m, labels, constraints, windows).
// Usage: ffbd_sim OUTFILE SEED REP [root|origin] -> writes a ffbd input (start = truth) and prints the truth line to stdout.
// Origin mode: one lineage from an origin age ~ Uniform(15, 25); fossil F2 has no clade constraint (it may sit anywhere, including above the root).
// Rate shifts: ffbd_sim OUTFILE SEED REP root|origin T1,T2,... simulates piecewise-constant rates (interval i = ages [T_{i-1}, T_i), i = 0
// youngest), with lambda_i, mu_i, psi_i ~ Exp(10) independently; the input gets rate_shifts and per-interval init values, and the truth
// line gets the per-interval rates appended (lambda_0..K, mu_0..K, psi_0..K). Without shifts the output is unchanged.
// Optional 6th argument, e.g. lambda,mu: only these rates vary across intervals; the others are one Exp(10) value for all intervals.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <random>
#include <string>
#include <vector>
using namespace std;
static mt19937_64 eng;
static double U() { return (eng() >> 11) * (1.0 / 9007199254740992.0); }
static double EXP(double rate) { double u; do { u = U(); } while (u <= 0); return -log(u) / rate; }
static double NORM(double mean, double sd) { double a; do { a = U(); } while (a <= 0); return mean + sd * sqrt(-2 * log(a)) * cos(6.283185307179586 * U()); }

struct N { int left = -1, right = -1; double age = 0; int kind = 0; }; // kind 0 extant tip, 1 fossil tip, 2 bifurcation, 3 SA node (right = continuing lineage, fossil tip created as left)
static vector<N> nodes; static int nExtant, nFossil;
static const int MAXE = 4, MAXF = 3;
static double lam, mu, psi_, rho = 0.5;
static vector<double> shiftT, lamV, muV, psiV; // piecewise-constant rates (shiftT empty: constant rates lam, mu, psi_)
static bool abortSim;

static int newNode(int kind, double age, int l = -1, int r = -1) { nodes.push_back({l, r, age, kind}); return (int)nodes.size() - 1; }
static int sim(double t) {
    if (abortSim) return -1;
    double a, l = lam, m = mu, p = psi_;
    if (shiftT.empty()) a = t - EXP(l + m + p);
    else { // waiting time across intervals (memoryless: on reaching a shift, draw again with the younger interval's rates)
        int i = 0; while (i < (int)shiftT.size() && t >= shiftT[i]) i++;
        a = t;
        while (true) {
            l = lamV[i]; m = muV[i]; p = psiV[i];
            double lo = i == 0 ? 0.0 : shiftT[i - 1], b = a - EXP(l + m + p);
            if (b > lo || i == 0) { a = b; break; }
            a = lo; i--;
        }
    }
    double rate = l + m + p;
    if (a <= 0) { if (U() < rho) { nExtant++; if (nExtant > MAXE) abortSim = true; return newNode(0, 0.0); } return -1; }
    double u = U() * rate;
    if (u < l) {
        int L = sim(a); int R = sim(a); if (L < 0) return R; if (R < 0) return L; return newNode(2, a, L, R);
    } else if (u < l + m) return -1;
    // fossil sample
    nFossil++; if (nFossil > MAXF) { abortSim = true; return -1; }
    int fid = newNode(1, a);
    int d = sim(a);
    if (d < 0) return fid; // fossil tip
    nodes[fid].kind = 1;
    return newNode(3, a, fid, d); // SA: left child is the fossil tip (zero-length), right child continues
}

int main(int argc, char **argv) {
    string out = argv[1]; eng.seed(stoull(argv[2])); int rep = argc > 3 ? atoi(argv[3]) : 0;
    bool originMode = (argc > 4 && string(argv[4]) == "origin");
    if (argc > 5) { string t; string a5 = argv[5]; size_t q = 0; while (q < a5.size()) { size_t e = a5.find(',', q); if (e == string::npos) e = a5.size(); shiftT.push_back(stod(a5.substr(q, e - q))); q = e + 1; } }
    int nI = (int)shiftT.size() + 1;
    bool varyW[3] = {true, true, true}; // optional 6th argument: which rates vary across intervals (as ffbd shift_rates); the others are shared
    if (argc > 6) { string a6 = argv[6]; varyW[0] = a6.find("lambda") != string::npos; varyW[1] = a6.find("mu") != string::npos; varyW[2] = a6.find("psi") != string::npos; }
    // fossil windows, label order fixed: F1 in [6,14], F2 in [2,16], F3 in [8,18]; labels: A,B,C,D,F1,F2,F3
    double wMin[3] = {6, 2, 8}, wMax[3] = {14, 16, 18};
    double se = 0.03;
    long tries = 0;
    while (true) {
        tries++;
        lam = EXP(10); mu = EXP(10); psi_ = EXP(10);
        if (!shiftT.empty()) { lamV.assign(nI, 0); muV.assign(nI, 0); psiV.assign(nI, 0); lamV[0] = lam; muV[0] = mu; psiV[0] = psi_; for (int i = 1; i < nI; i++) { lamV[i] = varyW[0] ? EXP(10) : lam; muV[i] = varyW[1] ? EXP(10) : mu; psiV[i] = varyW[2] ? EXP(10) : psi_; } }
        double rootAge = originMode ? 15.0 + 10.0 * U() : NORM(20, 2); if (rootAge < 1) continue;
        nodes.clear(); nExtant = nFossil = 0; abortSim = false;
        int rootId;
        if (originMode) {
            rootId = sim(rootAge);
            if (abortSim || rootId < 0 || nExtant != MAXE || nFossil != MAXF) continue;
            if (nodes[rootId].kind < 2) continue;
        } else {
            int rL = sim(rootAge); int rR = sim(rootAge);
            if (abortSim || rL < 0 || rR < 0 || nExtant != MAXE || nFossil != MAXF) continue;
            // build binary tree: node ids in `nodes`; root is a new node
            rootId = newNode(2, rootAge, rL, rR);
        }
        // gather tips and internal nodes (kind 3 = internal "SA parent" with fossil tip left child, zero branch length)
        vector<int> ext, fos;
        for (size_t i = 0; i < nodes.size(); i++) { if (nodes[i].kind == 0) ext.push_back((int)i); else if (nodes[i].kind == 1) fos.push_back((int)i); }
        if ((int)ext.size() != MAXE || (int)fos.size() != MAXF) continue;
        // parent links & masks
        int nn = (int)nodes.size(); vector<int> par(nn, -1);
        for (int i = 0; i < nn; i++) if (nodes[i].kind >= 2) { par[nodes[i].left] = i; par[nodes[i].right] = i; }
        // random labels: extant labels 0..3 (A,B,C,D), fossil labels 0..2 (F1,F2,F3); tip label index within mask: A=0..D=3, F1=4..F3=6
        vector<int> eperm = {0, 1, 2, 3}, fperm = {0, 1, 2};
        for (int i = 3; i > 0; i--) { int j = (int)(U() * (i + 1)); swap(eperm[i], eperm[j]); }
        for (int i = 2; i > 0; i--) { int j = (int)(U() * (i + 1)); swap(fperm[i], fperm[j]); }
        vector<int> tipLab(nn, -1);
        for (int i = 0; i < MAXE; i++) tipLab[ext[i]] = eperm[i];
        for (int i = 0; i < MAXF; i++) tipLab[fos[i]] = 4 + fperm[i];
        // windows
        bool ok = true; for (int i = 0; i < MAXF; i++) { double a = nodes[fos[i]].age; int w = fperm[i]; if (a < wMin[w] || a > wMax[w]) ok = false; }
        if (!ok) continue;
        // masks
        vector<int> mask(nn, 0);
        function<int(int)> rec = [&](int v) { if (nodes[v].kind < 2) { mask[v] = 1 << tipLab[v]; return mask[v]; } mask[v] = rec(nodes[v].left) | rec(nodes[v].right); return mask[v]; };
        rec(rootId);
        int emask = 0xF; bool ab = false, cd = false;
        int c1 = (1 << 0) | (1 << 1) | (1 << 4) | (1 << 6); // {A,B,F1,F3}
        int c2 = (1 << 2) | (1 << 3) | (1 << 5);            // {C,D,F2}
        bool hasC1 = false, hasC2 = false;
        for (int i = 0; i < nn; i++) { if (nodes[i].kind < 2) continue; if ((mask[i] & emask) == 0b0011) ab = true; if ((mask[i] & emask) == 0b1100) cd = true; if (mask[i] == c1) hasC1 = true; if (mask[i] == c2) hasC2 = true; }
        if (!ab || !cd || !hasC1 || (!originMode && !hasC2)) continue;
        // ACCEPT. Assign ffbd ids: tips A,B,C,D,F1,F2,F3 = 0..6; internal 7..12
        map<int, int> idOf; vector<int> labelToNode(7);
        for (int i = 0; i < nn; i++) if (nodes[i].kind < 2) labelToNode[tipLab[i]] = i;
        for (int l = 0; l < 7; l++) idOf[labelToNode[l]] = l;
        int nextId = 7; vector<int> internals;
        function<void(int)> order = [&](int v) { if (nodes[v].kind < 2) return; internals.push_back(v); idOf[v] = nextId++; order(nodes[v].left); order(nodes[v].right); };
        order(rootId);
        int nNodes = 13;
        // rates on the 6 extant branches: extant topology ((A,B)n1,(C,D)n2)n3 ; branches A,B,n1,C,D,n2 ; rates ~ Exp(mean m), m ~ Exp(10)
        double m = EXP(10); vector<double> rates(6); for (auto &r : rates) r = EXP(1.0 / m);
        // extant node ages: MRCA of extant labels pairs: node masks restricted to extant labels
        function<int(int)> extantMask = [&](int v) { if (nodes[v].kind == 0) return mask[v]; if (nodes[v].kind == 1) return 0; return extantMask(nodes[v].left) | extantMask(nodes[v].right); };
        // n1: extant mask {A,B} youngest node with both children containing extant; same for n2, n3 root of extant tree
        int n1 = -1, n2 = -1, n3 = -1;
        for (int v : internals) { int ml = extantMask(nodes[v].left), mr = extantMask(nodes[v].right); if (ml && mr) { int mm = ml | mr; if (mm == 0b0011) n1 = v; else if (mm == 0b1100) n2 = v; else if (mm == 0b1111) n3 = v; } }
        if (n1 < 0 || n2 < 0 || n3 < 0) { continue; }
        double tn1 = nodes[n1].age, tn2 = nodes[n2].age, tn3 = nodes[n3].age;
        // branches (parent ext id, child ext id or -1): A(n1,tip) B(n1,tip) n1(n3,n1) C(n2,tip) D(n2,tip) n2(n3,n2)
        vector<array<int, 2>> br = {{n1, -1}, {n1, -1}, {n3, n1}, {n2, -1}, {n2, -1}, {n3, n2}};
        auto T = [&](int b) { return nodes[br[b][0]].age - (br[b][1] < 0 ? 0.0 : nodes[br[b][1]].age); };
        // unrooted edges: the two root branches (n1, n2: branches 2 and 5) share one edge; terms: A,B,C,D,(n1+n2) = 5 terms
        vector<vector<int>> terms = {{0}, {1}, {3}, {4}, {2, 5}};
        vector<double> obs(5);
        for (int j = 0; j < 5; j++) { double mu_j = 0; for (int b : terms[j]) mu_j += rates[b] * T(b); obs[j] = NORM(fabs(mu_j), se); }
        // write input
        FILE *f = fopen(out.c_str(), "w");
        if (originMode) fprintf(f, "FFBD 1\ncond origin\nmode fast\nn 4\nm 3\nrho %.17g\nroot_prior none\norigin_prior uniform 15 25\n", rho);
        else fprintf(f, "FFBD 1\ncond root\nmode fast\nn 4\nm 3\nrho %.17g\nroot_prior normal 20 2\norigin_prior none\n", rho);
        if (shiftT.empty()) fprintf(f, "init_lambda %.17g\ninit_mu %.17g\ninit_psi %.17g\ninit_cond_age %.17g\ninit_clock %.17g\n", lam, mu, psi_, rootAge, m);
        else {
            fprintf(f, "rate_shifts"); for (double x : shiftT) fprintf(f, " %.17g", x); fprintf(f, "\n");
            if (argc > 6) fprintf(f, "shift_rates %s\n", argv[6]);
            const vector<double> *V[3] = {&lamV, &muV, &psiV}; const char *nm3[3] = {"init_lambda", "init_mu", "init_psi"};
            for (int w = 0; w < 3; w++) { fprintf(f, "%s", nm3[w]); for (double x : *V[w]) fprintf(f, " %.17g", x); fprintf(f, "\n"); }
            fprintf(f, "init_cond_age %.17g\ninit_clock %.17g\n", rootAge, m);
        }
        const char *nm[7] = {"A", "B", "C", "D", "FOSSIL_1", "FOSSIL_2", "FOSSIL_3"};
        fprintf(f, "tips 7\n");
        for (int i = 0; i < 4; i++) fprintf(f, "%d %s E 0 0\n", i, nm[i]);
        for (int i = 0; i < 3; i++) fprintf(f, "%d %s F %.17g %.17g\n", 4 + i, nm[4 + i], wMin[i], wMax[i]);
        fprintf(f, "initial_tree 13\n");
        // tips
        for (int l = 0; l < 7; l++) { int v = labelToNode[l]; int p = par[v]; int sa = (nodes[p].kind == 3 && nodes[v].kind == 1 && nodes[p].left == v) ? 1 : 0; fprintf(f, "%d %d -1 -1 %.17g %d\n", l, idOf[p], nodes[v].age, sa); }
        for (int v : internals) { int p = v == rootId ? -1 : par[v]; int pid = p < 0 ? -1 : idOf[p]; fprintf(f, "%d %d %d %d %.17g 0\n", idOf[v], pid, idOf[nodes[v].left], idOf[nodes[v].right], nodes[v].age); }
        if (originMode) fprintf(f, "constraints 1\ncl_ABF13 4 0 1 4 6\n"); else fprintf(f, "constraints 2\ncl_ABF13 4 0 1 4 6\ncl_CDF2 3 2 3 5\n");
        fprintf(f, "extant_nodes 3\n%d n1\n%d n2\n%d n3\n", idOf[n1], idOf[n2], idOf[n3]);
        fprintf(f, "branches 6\n");
        for (int b = 0; b < 6; b++) fprintf(f, "%d %d %d\n", b, idOf[br[b][0]], br[b][1] < 0 ? -1 : idOf[br[b][1]]);
        fprintf(f, "terms 5\n");
        for (int j = 0; j < 5; j++) { fprintf(f, "%d %.17g %.17g %zu", j, obs[j], se, terms[j].size()); for (int b : terms[j]) fprintf(f, " %d", b); fprintf(f, "\n"); }
        fclose(f);
        // rates are start = truth; the sampler starts rates at m, which is fine for burn-in; write them for the truth record
        double meanR = 0; for (double r : rates) meanR += r / 6.0;
        printf("TRUTH\t%d\t%ld\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g", rep, tries, lam, mu, psi_, rootAge, m, tn1, tn2, tn3, nodes[labelToNode[4]].age, nodes[labelToNode[5]].age, nodes[labelToNode[6]].age, meanR);
        int nsa = 0; for (int i = 0; i < nn; i++) if (nodes[i].kind == 3) nsa++;
        printf("\t%d\t%.10g", nsa, nodes[rootId].age);
        if (!shiftT.empty()) { for (double x : lamV) printf("\t%.10g", x); for (double x : muV) printf("\t%.10g", x); for (double x : psiV) printf("\t%.10g", x); }
        printf("\n");
        return 0;
    }
}
