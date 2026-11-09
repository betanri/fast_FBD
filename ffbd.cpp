// ffbd.cpp -- dedicated MCMC sampler for the Fast FBD model of FBD/tree/fast_fbd_test/rev/fbd_run_root.Rev.
// C++17, standard library only. See README.md. Input file written by scripts/ffbd_convert.py.
//
// Model (matches RevBayes 1.3.2 dnBDSTP(rootAge|originAge, lambda, mu, psi, rho, r=0, taxa) inside dnConstrainedTopology):
// - FBD density as in BirthDeathSamplingTreatmentProcess::computeLnProbabilityTimes + lnProbTreeShape.
// - Each fossil age is uniform in its window (density 1/(max-min+2e-4), the same constant as the Rev window terms).
// - Mode "fast": rate_i ~ Exp(mean m), m ~ Exp(10), b_obs_j ~ Normal(|sum rate_i*time_i|, se_j) (edges of the unrooted tree).
// - Mode "control": no rates, no data.
#include <algorithm>
#include <array>
#include <chrono>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
static const double LN2 = 0.6931471805599453;
static const double LN2PI = 1.8378770664093453;

struct RNG {
    mt19937_64 e;
    explicit RNG(uint64_t s = 1) : e(s) {}
    inline double u() { return (e() >> 11) * (1.0 / 9007199254740992.0); }
    inline int below(int n) { return (int)(u() * n); }
    inline double normal() {
        double a = u(), b = u();
        if (a <= 0) a = 1e-300;
        return sqrt(-2 * log(a)) * cos(6.283185307179586 * b);
    }
};

// ---- FBD quantities (constant rate, r = 0), same formulae as RevBayes ----
struct FBD {
    double lam = 0, mu = 0, psi = 0, rho = 0, A = 0, B = 0, lnLam = 0, lnPsi = 0;
    void set(double l, double m, double p, double r) {
        lam = l; mu = m; psi = p; rho = r;
        A = sqrt((l - m - p) * (l - m - p) + 4 * l * p);
        B = ((1.0 - 2.0 * (1.0 - r)) * l + m + p) / A;
        lnLam = log(l); lnPsi = log(p);
    }
    inline double lnD(double t) const {
        if (t < DBL_EPSILON) return 0.0;
        double e = exp(-A * t);
        return 2.0 * LN2 - A * t - 2.0 * log(1.0 + B + e * (1.0 - B));
    }
    inline double lnE(double t) const {
        double e = exp(-A * t);
        double r = (1.0 + B - e * (1.0 - B)) / (1.0 + B + e * (1.0 - B));
        return log((lam + mu + psi - A * r) / (2.0 * lam));
    }
};

static inline double sq(double x) { return x * x; }

// ---- root prior ----
struct RootPrior {
    int type = 0; // 0 none (uniform improper never used), 1 normal(mean, sd), 2 uniform(min,max), 3 soft(min,max,sd,p)
    double a = 0, b = 0, sd = 0, p = 0;
    double ln(double x) const {
        if (type == 1) return -0.5 * LN2PI - log(sd) - 0.5 * sq((x - a) / sd);
        if (type == 2) return (x < a || x > b) ? -INFINITY : -log(b - a);
        if (type == 3) {
            // RevBayes SoftBoundUniformNormalDistribution::computeLnProbability, unchanged (note its log(1 - p/2))
            if (x < a) return log(1.0 - p / 2.0) + (-0.5 * LN2PI - log(sd) - 0.5 * sq((x - a) / sd));
            if (x > b) return log(1.0 - p / 2.0) + (-0.5 * LN2PI - log(sd) - 0.5 * sq((x - b) / sd));
            return log(p) - log(b - a);
        }
        return 0.0;
    }
};

// ---- tunable proposal parameter ----
struct Tune {
    double val = 0.3; int att = 0, acc = 0; long totalAtt = 0, totalAcc = 0; int cycles = 0;
    void hit(bool a) { att++; if (a) acc++; totalAtt++; if (a) totalAcc++; }
    void adapt(double target, double lo, double hi) {
        if (att < 100) return;
        double r = (double)acc / att; cycles++;
        double step = (r - target) / sqrt((double)cycles);
        val *= exp(step); val = min(hi, max(lo, val));
        att = 0; acc = 0;
    }
};

enum MoveType { M_NODE, M_COUPLED, M_FOSAGE, M_REGRAFT_L, M_REGRAFT_G, M_COLLEXP, M_ROOT, M_ROOTC, M_ORIGIN,
                M_RATE, M_RATESCALE_JOINT, M_RATESCALE_ALL, M_MEAN, M_SIGMA, M_LAMBDA, M_MU, M_PSI, M_UPDOWN, N_MOVES };
static const char *MOVE_NAMES[] = {"node_slide", "node_slide_coupled", "fossil_age", "regraft_local", "regraft_global",
    "collapse_expand", "root_slide", "root_slide_coupled", "origin_slide", "rate_scale", "rates_mean_joint_scale",
    "rates_vector_scale", "mean_rate_scale", "ucln_sigma_scale", "lambda_scale", "mu_scale", "psi_scale", "lambda_mu_updown"};

struct Options {
    string input, prefix = "ffbd_run";
    uint64_t seed = 1; long gens = 1000000, printEvery = 1000, treeEvery = 10000, ckpEvery = 0, burnin = 0;
    string mode = ""; // "", "fast", "control"
    bool resume = false, check = false, noData = false, rbSA = false, logRates = false, logFossils = false, fixFbd = false; bool checkPerStep = false;
    string dumpFile = ""; long dumpEvery = 0; string noMoves = "";
};

class Chain {
public:
    int n = 0, mF = 0, N = 0, nNodes = 0, root = -1;
    bool useOrigin = false; string cond, mode; bool hasClock = true, useData = true;
    double rho = 1.0, lnRho = 0.0; RootPrior rootPrior; RootPrior originPrior;
    double origin = 0.0; // origin age (origin mode); unused in root mode

    // tree arrays
    vector<int> par, c0, c1; vector<double> age; vector<char> sa; vector<string> names;
    vector<double> winMin, winMax; double lnWin = 0.0;

    // extant bookkeeping
    vector<char> isExtInt; vector<int> extNode; vector<string> extName; // extant internal nodes
    double uclnMu = 0; int K = 0, J = 0;                                                   // branches, terms
    vector<int> bpar, bch; vector<vector<int>> tb; vector<double> obs, se, lnSe; vector<array<int, 4>> nodeTerms; vector<int> nNodeTerms;
    vector<array<int, 3>> nodeBranches; vector<int> nNodeBranches;
    vector<double> tll; vector<int> bTerm;

    // constraints
    int nCons = 0; vector<vector<int>> consMembers; vector<int> parentCons; vector<int> x; vector<int> homeC; vector<vector<int>> region; vector<int> label, pos;
    vector<string> consName;

    // parameters
    double lam = 0, mu = 0, psi = 0, m = 0; vector<double> rates; double sumRates = 0;
    bool ucln = false; double sigma = 0.5; vector<double> calLo, calHi; // UCLN clock (catfish) and uniform calibrations on extant node ages
    FBD P;

    // cached log densities
    double lnFBD = 0, lnRoot = 0, lnHyper = 0, lnRatePrior = 0, lnData = 0;
    double lnCalConst = 0.0; vector<double> initRates;
    long iter = 0; long burnDone = 0; RNG rng;
    Tune tLam, tMu, tPsi, tUD, tRootS, tRootCo, tOrigin, tRate, tMean, tJoint, tVec, tSigma;
    array<long, N_MOVES> att{}, acc{};
    vector<double> weights; vector<double> cum; double totalW = 0;
    Options opt;
    bool tuning = false;

    explicit Chain(const Options &o) : rng(o.seed), opt(o) {}

    void read(const string &path);
    void init(); void setupMoves();
    double lnFBDFull(const FBD &p) const; void recomputeAll();
    double lnHyperFull() const { return (log(10.0) - 10.0 * lam) + (log(10.0) - 10.0 * mu) + (log(10.0) - 10.0 * psi); }
    inline double top(int v) const { return par[v] == -1 ? origin : age[par[v]]; }
    inline int sib(int v) const { int p = par[v]; return c0[p] == v ? c1[p] : c0[p]; }
    inline bool isSAParent(int v) const { return v >= N && ((c0[v] < N && sa[c0[v]]) || (c1[v] < N && sa[c1[v]])); }
    inline double branchTime(int i) const { double pa = age[bpar[i]], ca = bch[i] < 0 ? 0.0 : age[bch[i]]; return fabs(pa - ca); }
    double termMean(int j) const { double s = 0; for (int b : tb[j]) s += rates[b] * branchTime(b); return fabs(s); }
    inline double termLL(int j) const { return -0.5 * LN2PI - lnSe[j] - 0.5 * sq((obs[j] - termMean(j)) / se[j]); }
    inline double lnRateP1(double r) const {
        if (ucln) { double z = (log(r) - uclnMu) / sigma; return -log(r) - log(sigma) - 0.5 * LN2PI - 0.5 * z * z; }
        return log(1.0 / m) - r / m;
    }
    double rateLnPrior() const { double s = 0; for (double r : rates) s += lnRateP1(r); s += log(10.0) - 10.0 * m; if (ucln) s += log(3.0) - 3.0 * sigma; return s; }
    double rateLnPriorGiven(const vector<double> &rs, double mNew, double sNew) const {
        double uMu = log(mNew) - 0.5 * sNew * sNew, s = 0;
        for (double r : rs) { if (ucln) { double z = (log(r) - uMu) / sNew; s += -log(r) - log(sNew) - 0.5 * LN2PI - 0.5 * z * z; } else s += log(1.0 / mNew) - r / mNew; }
        s += log(10.0) - 10.0 * mNew; if (ucln) s += log(3.0) - 3.0 * sNew; return s;
    }
    void setLabel(int v, int L);
    void prune(int p, int s, int g); void attach(int p, int f, int c, double a);
    int pickRegion(int C, int p, int f);
    bool eligible(int x, int C, int p, int f) const;
    int neighbors(int c, int C, int p, int f, array<int, 4> &out) const;

    // moves
    bool moveNode(); bool moveCoupled(); bool moveFossilAge(); bool moveRegraft(bool local); bool moveCollapseExpand();
    bool moveRoot(bool coupled); bool moveOrigin(); bool moveRate(); bool moveJoint(); bool moveVector(); bool moveMean(); bool moveSigma();
    bool moveFBDScale(int which); bool moveUpDown();
    void validate();

    // output
    void run(); void writeHeader(FILE *lf, FILE *tf); void writeLog(FILE *lf); void writeTree(FILE *tf);
    string newick() const; void saveCheckpoint(const string &p) const; void loadCheckpoint(const string &p);
    void dumpState(FILE *f);
    void summary() const;
};

// ------------------------------------------------------------------------------------------------
static void die(const string &m) { cerr << "ERROR: " << m << endl; exit(2); }

void Chain::read(const string &path) {
    ifstream in(path); if (!in) die("cannot open input " + path);
    string tok; int fmt; in >> tok >> fmt; if (tok != "FFBD") die("not a ffbd input file");
    map<string, string> kv;
    string key;
    // scalar header up to the "tips" keyword
    string rootp, origp; vector<double> rootv, origv;
    double initLam = 0, initMu = 0, initPsi = 0, initCond = 0, initClock = 0;
    while (in >> key) {
        if (key == "tips") break;
        if (key == "cond") in >> cond;
        else if (key == "mode") { in >> mode; }
        else if (key == "n") in >> n;
        else if (key == "m") in >> mF;
        else if (key == "rho") in >> rho;
        else if (key == "root_prior") { in >> rootp; int c = rootp == "normal" ? 2 : rootp == "uniform" ? 2 : rootp == "soft" ? 4 : 0; rootv.resize(c); for (auto &v : rootv) in >> v; }
        else if (key == "origin_prior") { in >> origp; int c = origp == "uniform" ? 2 : 0; origv.resize(c); for (auto &v : origv) in >> v; }
        else if (key == "init_lambda") in >> initLam;
        else if (key == "init_mu") in >> initMu;
        else if (key == "init_psi") in >> initPsi;
        else if (key == "init_cond_age") in >> initCond;
        else if (key == "init_clock") in >> initClock;
        else if (key == "rate_prior") { string rp; in >> rp; ucln = (rp == "ucln"); }
        else if (key == "init_sigma") in >> sigma;
        else die("unknown header key " + key);
    }
    useOrigin = (cond == "origin");
    if (rootp == "normal") { rootPrior.type = 1; rootPrior.a = rootv[0]; rootPrior.sd = rootv[1]; }
    else if (rootp == "uniform") { rootPrior.type = 2; rootPrior.a = rootv[0]; rootPrior.b = rootv[1]; }
    else if (rootp == "soft") { rootPrior.type = 3; rootPrior.a = rootv[0]; rootPrior.b = rootv[1]; rootPrior.sd = rootv[2]; rootPrior.p = rootv[3]; }
    if (origp == "uniform") { originPrior.type = 2; originPrior.a = origv[0]; originPrior.b = origv[1]; }
    N = n + mF; nNodes = 2 * N - 1;
    int cnt; in >> cnt; if (cnt != N) die("tips count");
    names.resize(N); winMin.resize(N); winMax.resize(N);
    for (int i = 0; i < N; i++) { int id; string nm; string typ; double mn, mx; in >> id >> nm >> typ >> mn >> mx; if (id != i) die("tip order"); names[i] = nm; winMin[i] = mn; winMax[i] = mx; }
    in >> key >> cnt; if (key != "initial_tree" || cnt != nNodes) die("initial_tree section " + key);
    par.assign(nNodes, -1); c0.assign(nNodes, -1); c1.assign(nNodes, -1); age.assign(nNodes, 0.0); sa.assign(nNodes, 0);
    for (int i = 0; i < nNodes; i++) { int id, p, a, b; double t; int s; in >> id >> p >> a >> b >> t >> s; par[id] = p; c0[id] = a; c1[id] = b; age[id] = t; sa[id] = (char)s; if (p == -1) root = id; }
    in >> key >> nCons; if (key != "constraints") die("constraints section");
    // constraint 0 is the pseudo-constraint "whole tree"
    consMembers.assign(nCons + 1, {}); consName.assign(nCons + 1, "ROOT"); consMembers[0].resize(N); iota(consMembers[0].begin(), consMembers[0].end(), 0);
    for (int k = 1; k <= nCons; k++) { int nm; string name; in >> name >> nm; consName[k] = name; consMembers[k].resize(nm); for (auto &v : consMembers[k]) in >> v; }
    in >> key >> cnt; if (key != "extant_nodes") die("extant_nodes section");
    extNode.resize(cnt); extName.resize(cnt);
    for (int i = 0; i < cnt; i++) { in >> extNode[i] >> extName[i]; }
    in >> key >> K; if (key != "branches") die("branches section");
    bpar.resize(K); bch.resize(K);
    for (int i = 0; i < K; i++) { int id, a, b; in >> id >> a >> b; bpar[id] = a; bch[id] = b; }
    in >> key >> J; if (key != "terms") die("terms section");
    tb.resize(J); obs.resize(J); se.resize(J); lnSe.resize(J);
    for (int j = 0; j < J; j++) { int id, nb; double o, s; in >> id >> o >> s >> nb; obs[id] = o; se[id] = s; lnSe[id] = log(s); tb[id].resize(nb); for (auto &v : tb[id]) in >> v; }
    vector<double> initR;
    calLo.assign(nNodes, -1.0); calHi.assign(nNodes, -1.0);
    while (in >> key) {
        if (key == "calibrations") { int c; in >> c; for (int i = 0; i < c; i++) { int v; double lo, hi; in >> v >> lo >> hi; calLo[v] = lo; calHi[v] = hi; lnCalConst -= log(hi - lo); } }
        else if (key == "init_rates") { int c; in >> c; initR.resize(c); for (auto &v : initR) in >> v; }
        else die("unknown section after terms: " + key);
    }
    initRates = initR;
    lam = initLam; mu = initMu; psi = initPsi; m = initClock; if (useOrigin) origin = initCond; else age[root] = initCond;
    if (!opt.mode.empty()) mode = opt.mode;
    hasClock = (mode == "fast"); useData = hasClock && !opt.noData;
    lnRho = log(rho);
}

void Chain::init() {
    if (!hasClock) { K = 0; J = 0; bpar.clear(); bch.clear(); tb.clear(); obs.clear(); se.clear(); lnSe.clear(); }
    rates.assign(K, m); if ((int)initRates.size() == K) rates = initRates; sumRates = accumulate(rates.begin(), rates.end(), 0.0);
    uclnMu = log(m) - 0.5 * sigma * sigma;
    // fossil windows: constant window density (Rev: 1/(max-min+2e-4)), hard bounds in the model
    lnWin = lnCalConst; for (int f = n; f < N; f++) lnWin -= log(winMax[f] - winMin[f] + 2e-4);
    // extant internal nodes and clock incidence
    isExtInt.assign(nNodes, 0); for (int v : extNode) isExtInt[v] = 1;
    nodeTerms.assign(nNodes, array<int, 4>{-1, -1, -1, -1}); nNodeTerms.assign(nNodes, 0);
    nodeBranches.assign(nNodes, array<int, 3>{-1, -1, -1}); nNodeBranches.assign(nNodes, 0);
    vector<int> termOfBranch(K, -1); for (int j = 0; j < J; j++) for (int b : tb[j]) termOfBranch[b] = j;
    bTerm = termOfBranch;
    for (int b = 0; b < K; b++) for (int v : {bpar[b], bch[b]}) if (v >= 0) {
        if (nNodeBranches[v] >= 3) die("extant node with >3 incident branches");
        nodeBranches[v][nNodeBranches[v]++] = b;
    }
    for (int v = 0; v < nNodes; v++) for (int i = 0; i < nNodeBranches[v]; i++) {
        int j = termOfBranch[nodeBranches[v][i]]; bool seen = false; for (int q = 0; q < nNodeTerms[v]; q++) if (nodeTerms[v][q] == j) seen = true;
        if (!seen) nodeTerms[v][nNodeTerms[v]++] = j;
    }
    tll.assign(J, 0.0);
    // constraints: laminar check, parent constraint, home constraint of each fossil, node of each constraint, labels
    vector<vector<char>> mem(nCons + 1, vector<char>(N, 0));
    for (int k = 0; k <= nCons; k++) for (int v : consMembers[k]) mem[k][v] = 1;
    parentCons.assign(nCons + 1, 0); parentCons[0] = -1;
    for (int k = 1; k <= nCons; k++) {
        int best = 0; size_t bs = N + 1;
        for (int q = 1; q <= nCons; q++) if (q != k && consMembers[q].size() > consMembers[k].size()) {
            bool sub = true; for (int v : consMembers[k]) if (!mem[q][v]) { sub = false; break; }
            if (sub && consMembers[q].size() < bs) { best = q; bs = consMembers[q].size(); }
        }
        parentCons[k] = best;
    }
    homeC.assign(N, 0);
    for (int f = n; f < N; f++) {
        int best = 0; size_t bs = N + 1; for (int k = 1; k <= nCons; k++) if (mem[k][f] && consMembers[k].size() < bs) { best = k; bs = consMembers[k].size(); }
        homeC[f] = best;
        if (best == 0 && !useOrigin) die("root-conditioned model: fossil " + names[f] + " has no constraint below the root (not supported, see README)");
    }
    // subtree tips for labelling
    vector<vector<int>> kids(nNodes); vector<int> order; order.reserve(nNodes);
    for (int v = 0; v < nNodes; v++) if (v >= N) kids[v] = {c0[v], c1[v]};
    vector<int> st{root};
    while (!st.empty()) { int v = st.back(); st.pop_back(); order.push_back(v); for (int k : kids[v]) st.push_back(k); }
    if ((int)order.size() != nNodes) die("initial tree is not connected");
    vector<vector<int>> tips(nNodes);
    for (auto it = order.rbegin(); it != order.rend(); ++it) { int v = *it; if (v < N) tips[v] = {v}; else { tips[v] = tips[c0[v]]; tips[v].insert(tips[v].end(), tips[c1[v]].begin(), tips[c1[v]].end()); sort(tips[v].begin(), tips[v].end()); } }
    x.assign(nCons + 1, -1);
    for (int k = 1; k <= nCons; k++) {
        int found = -1; vector<int> s = consMembers[k]; sort(s.begin(), s.end());
        for (int v = 0; v < nNodes; v++) if (tips[v].size() == s.size() && tips[v] == s) { found = v; break; }
        if (found < 0) die("initial tree violates constraint " + consName[k]);
        x[k] = found;
    }
    region.assign(nCons + 1, {}); label.assign(nNodes, -1); pos.assign(nNodes, -1);
    for (int v = 0; v < nNodes; v++) {
        if (v == root && !useOrigin) continue; // no edge above the root
        int L = 0; size_t bs = N + 1;
        for (int k = 1; k <= nCons; k++) if (consMembers[k].size() < bs && consMembers[k].size() > tips[v].size()) {
            bool sub = true; for (int q : tips[v]) if (!mem[k][q]) { sub = false; break; }
            if (sub) { L = k; bs = consMembers[k].size(); }
        }
        for (int k = 1; k <= nCons; k++) if (x[k] == v) L = parentCons[k]; // the node of constraint k sits in its parent's region
        setLabel(v, L);
    }
    setupMoves();
    recomputeAll();
    tLam.val = tMu.val = tPsi.val = 0.3; tUD.val = 0.2; tRootS.val = 1.0; tRootCo.val = 1.0; tOrigin.val = 2.0; tRate.val = 0.5; tMean.val = 0.3; tJoint.val = 0.1; tVec.val = 0.05; tSigma.val = 0.2;
}

void Chain::setLabel(int v, int L) {
    if (label[v] == L) return;
    if (label[v] >= 0) { auto &r = region[label[v]]; int i = pos[v]; int b = r.back(); r[i] = b; pos[b] = i; r.pop_back(); }
    label[v] = L;
    if (L >= 0) { auto &r = region[L]; pos[v] = (int)r.size(); r.push_back(v); } else pos[v] = -1;
}

double Chain::lnFBDFull(const FBD &p) const {
    double s = 0; int nsa = 0;
    for (int f = n; f < N; f++) {
        s -= p.lnD(age[f]);
        if (sa[f]) { s += p.lnPsi; nsa++; } else s += p.lnPsi + p.lnE(age[f]);
    }
    for (int v = N; v < nNodes; v++) {
        s += p.lnD(age[v]);
        if (!isSAParent(v) && (useOrigin || v != root)) s += p.lnLam;
    }
    if (useOrigin) s += p.lnD(origin); else s += p.lnD(age[root]);
    s += n * log(p.rho);
    s += (N - nsa - 1) * LN2 - lgamma(N + 1.0);
    return s;
}

void Chain::recomputeAll() {
    P.set(lam, mu, psi, rho);
    lnFBD = lnFBDFull(P);
    lnRoot = useOrigin ? originPrior.ln(origin) : rootPrior.ln(age[root]);
    lnHyper = lnHyperFull();
    sumRates = accumulate(rates.begin(), rates.end(), 0.0);
    lnRatePrior = hasClock ? rateLnPrior() : 0.0;
    lnData = 0; for (int j = 0; j < J; j++) { tll[j] = useData ? termLL(j) : 0.0; lnData += tll[j]; }
}

void Chain::prune(int p, int s, int g) {
    par[s] = g;
    if (g == -1) root = s; else { if (c0[g] == p) c0[g] = s; else c1[g] = s; }
}
void Chain::attach(int p, int f, int c, double a) {
    int g2 = par[c]; par[c] = p; par[f] = p; c0[p] = f; c1[p] = c; par[p] = g2; age[p] = a;
    if (g2 == -1) root = p; else { if (c0[g2] == c) c0[g2] = p; else c1[g2] = p; }
}

// ---- proposals --------------------------------------------------------------------------------
bool Chain::moveNode() {
    int v = N + rng.below(N - 1);
    if (v == root && !useOrigin) return false;
    if (isSAParent(v)) return false;
    double lo = max(age[c0[v]], age[c1[v]]), hi = top(v);
    if (calLo[v] >= 0) { lo = max(lo, calLo[v]); hi = min(hi, calHi[v]); }
    if (hi <= lo) return false;
    double oldA = age[v], newA = lo + rng.u() * (hi - lo);
    double dF = P.lnD(newA) - P.lnD(oldA);
    double dD = 0;
    array<double, 4> nv{};
    age[v] = newA;
    if (useData && nNodeTerms[v] > 0) for (int i = 0; i < nNodeTerms[v]; i++) { int j = nodeTerms[v][i]; nv[i] = termLL(j); dD += nv[i] - tll[j]; }
    if (log(rng.u()) < dF + dD) {
        lnFBD += dF; lnData += dD;
        if (useData) for (int i = 0; i < nNodeTerms[v]; i++) tll[nodeTerms[v][i]] = nv[i];
        return true;
    }
    age[v] = oldA; return false;
}

bool Chain::moveCoupled() {
    if (!hasClock || K == 0) return false;
    int v = extNode[rng.below((int)extNode.size())];
    if (v == root && !useOrigin) return false;
    if (isSAParent(v)) return false;
    double lo = max(age[c0[v]], age[c1[v]]), hi = top(v);
    if (calLo[v] >= 0) { lo = max(lo, calLo[v]); hi = min(hi, calHi[v]); }
    if (hi <= lo) return false;
    double oldA = age[v], newA = lo + rng.u() * (hi - lo);
    double dF = P.lnD(newA) - P.lnD(oldA);
    int nb = nNodeBranches[v]; array<double, 3> oldT{}, newT{}, oldR{}, newR{};
    for (int i = 0; i < nb; i++) { oldT[i] = branchTime(nodeBranches[v][i]); }
    age[v] = newA;
    double hast = 0, dP = 0;
    for (int i = 0; i < nb; i++) {
        int b = nodeBranches[v][i]; newT[i] = branchTime(b);
        if (newT[i] <= 0 || oldT[i] <= 0) { age[v] = oldA; return false; }
        oldR[i] = rates[b]; newR[i] = oldR[i] * oldT[i] / newT[i];
        hast += log(oldT[i] / newT[i]); dP += lnRateP1(newR[i]) - lnRateP1(oldR[i]);
    }
    if (log(rng.u()) < dF + dP + hast) {
        lnFBD += dF; lnRatePrior += dP;
        for (int i = 0; i < nb; i++) { rates[nodeBranches[v][i]] = newR[i]; sumRates += newR[i] - oldR[i]; }
        return true;
    }
    age[v] = oldA; return false;
}

bool Chain::moveFossilAge() {
    if (mF == 0) return false;
    int f = n + rng.below(mF);
    double lo = winMin[f], hi = winMax[f];
    if (hi <= lo) return false;
    if (sa[f]) {
        int p = par[f]; int s = sib(f);
        lo = max(lo, age[s]); hi = min(hi, top(p));
        if (hi <= lo) return false;
        double a = lo + rng.u() * (hi - lo); age[f] = a; age[p] = a; return true;
    }
    hi = min(hi, age[par[f]]); if (hi <= lo) return false;
    double oldA = age[f], newA = lo + rng.u() * (hi - lo);
    double dF = (P.lnE(newA) - P.lnD(newA)) - (P.lnE(oldA) - P.lnD(oldA));
    if (log(rng.u()) < dF) { age[f] = newA; lnFBD += dF; return true; }
    return false;
}

int Chain::pickRegion(int C, int p, int f) {
    const auto &arr = region[C];
    bool hasF = (label[f] == C), hasP = (label[p] == C);
    int k = (int)arr.size() - (hasF ? 1 : 0) - (hasP ? 1 : 0);
    int extra = (C != 0 && x[C] != p) ? 1 : 0;
    int total = k + extra;
    if (total <= 0) return -1;
    int r = rng.below(total);
    if (r == k) return x[C];
    int c; do { c = arr[rng.below((int)arr.size())]; } while (c == f || c == p);
    return c;
}

bool Chain::eligible(int y, int C, int p, int f) const {
    if (y == p || y == f) return false;
    if (y == root && !useOrigin) return false;
    if (C != 0 && y == x[C]) return true;
    return label[y] == C;
}

int Chain::neighbors(int c, int C, int p, int f, array<int, 4> &out) const {
    int k = 0;
    int pc = par[c];
    if (pc != -1) {
        // edge above pc (the root has one only in origin mode)
        if ((par[pc] != -1 || useOrigin) && eligible(pc, C, p, f)) out[k++] = pc;
        int sb = (c0[pc] == c) ? c1[pc] : c0[pc];
        if (eligible(sb, C, p, f)) out[k++] = sb;
    }
    if (c >= N) { if (eligible(c0[c], C, p, f)) out[k++] = c0[c]; if (eligible(c1[c], C, p, f)) out[k++] = c1[c]; }
    return k;
}

bool Chain::moveRegraft(bool local) {
    if (mF == 0) return false;
    int f = n + rng.below(mF);
    int p = par[f], s = (c0[p] == f) ? c1[p] : c0[p], g = par[p];
    if (p == root && !useOrigin) return false;
    if (s < N && sa[s]) return false;
    int C = homeC[f]; bool saF = sa[f] != 0; double ageF = age[f], oldAgeP = age[p];
    int A = (C == 0) ? -1 : (x[C] == p ? s : x[C]);
    int c = -1; double hast = 0.0;
    if (!local) { c = pickRegion(C, p, f); if (c < 0) return false; }
    prune(p, s, g);
    if (local) {
        array<int, 4> e1{}, e2{}; int ns = neighbors(s, C, p, f, e1);
        if (ns == 0) { attach(p, f, s, oldAgeP); return false; }
        c = e1[rng.below(ns)];
        int nc = neighbors(c, C, p, f, e2);
        hast += log((double)ns) - log((double)nc);
    }
    double topC = par[c] == -1 ? origin : age[par[c]];
    double dF = 0, newA = oldAgeP, logR = 0;
    if (saF) {
        if (!(age[c] < ageF && ageF < topC)) { attach(p, f, s, oldAgeP); return false; }
        newA = ageF; logR = hast;
    } else {
        double lo = max(age[c], ageF), hi = topC;
        if (hi <= lo) { attach(p, f, s, oldAgeP); return false; }
        newA = lo + rng.u() * (hi - lo);
        double rangeOld = (par[s] == -1 ? origin : age[par[s]]) - max(age[s], ageF);
        dF = P.lnD(newA) - P.lnD(oldAgeP);
        logR = dF + hast + log(hi - lo) - log(rangeOld);
    }
    if (log(rng.u()) < logR) {
        attach(p, f, c, newA);
        lnFBD += dF;
        if (C != 0) {
            int xNew = (c == A) ? p : A;
            x[C] = xNew;
            setLabel(p, xNew == p ? parentCons[C] : C);
            setLabel(A, xNew == A ? parentCons[C] : C);
        }
        return true;
    }
    attach(p, f, s, oldAgeP);
    return false;
}

bool Chain::moveCollapseExpand() {
    if (mF == 0) return false;
    int f = n + rng.below(mF);
    int p = par[f]; int s = sib(f); double t = age[f];
    if (sa[f]) {
        double hi = top(p), lo = opt.rbSA ? winMax[f] : t;
        if (hi <= lo) return false;
        double a = lo + rng.u() * (hi - lo);
        double dF = P.lnLam + P.lnE(t) + P.lnD(a) - P.lnD(t) + LN2;
        if (log(rng.u()) < dF + log(hi - lo)) { age[p] = a; sa[f] = 0; lnFBD += dF; return true; }
        return false;
    }
    if (p == root && !useOrigin) return false;
    if (age[s] >= t) return false;
    if (s < N && sa[s]) return false;
    double hi = top(p), lo = opt.rbSA ? winMax[f] : t;
    if (hi <= lo) return false;
    double dF = -(P.lnLam + P.lnE(t) + P.lnD(age[p]) - P.lnD(t) + LN2);
    if (log(rng.u()) < dF - log(hi - lo)) { age[p] = t; sa[f] = 1; lnFBD += dF; return true; }
    return false;
}

bool Chain::moveRoot(bool coupled) {
    if (useOrigin) return false; // origin mode: the root is an ordinary node
    Tune &tn = coupled ? tRootCo : tRootS;
    double oldA = age[root], newA = oldA + (rng.u() - 0.5) * 2.0 * tn.val;
    double lo = max(age[c0[root]], age[c1[root]]);
    if (newA <= lo) return false;
    if (calLo[root] >= 0 && (newA < calLo[root] || newA > calHi[root])) return false;
    double dP0 = rootPrior.ln(newA) - lnRoot;
    if (!isfinite(dP0)) return false;
    double dF = 2.0 * (P.lnD(newA) - P.lnD(oldA));
    bool ext = isExtInt[root] != 0;
    if (coupled && ext && hasClock) {
        int nb = nNodeBranches[root]; array<double, 3> oT{}, nT{}, oR{}, nR{};
        for (int i = 0; i < nb; i++) oT[i] = branchTime(nodeBranches[root][i]);
        age[root] = newA;
        double hast = 0, dP = 0;
        for (int i = 0; i < nb; i++) { int b = nodeBranches[root][i]; nT[i] = branchTime(b); oR[i] = rates[b]; nR[i] = oR[i] * oT[i] / nT[i]; hast += log(oT[i] / nT[i]); dP += lnRateP1(nR[i]) - lnRateP1(oR[i]); }
        if (log(rng.u()) < dF + dP0 + dP + hast) {
            lnFBD += dF; lnRoot += dP0; lnRatePrior += dP;
            for (int i = 0; i < nb; i++) { rates[nodeBranches[root][i]] = nR[i]; sumRates += nR[i] - oR[i]; }
            return true;
        }
        age[root] = oldA; return false;
    }
    age[root] = newA;
    double dD = 0; array<double, 4> nv{};
    if (useData && ext) for (int i = 0; i < nNodeTerms[root]; i++) { int j = nodeTerms[root][i]; nv[i] = termLL(j); dD += nv[i] - tll[j]; }
    if (log(rng.u()) < dF + dP0 + dD) {
        lnFBD += dF; lnRoot += dP0; lnData += dD;
        if (useData && ext) for (int i = 0; i < nNodeTerms[root]; i++) tll[nodeTerms[root][i]] = nv[i];
        return true;
    }
    age[root] = oldA; return false;
}

bool Chain::moveOrigin() {
    if (!useOrigin) return false;
    double oldO = origin, newO = oldO + (rng.u() - 0.5) * 2.0 * tOrigin.val;
    if (newO <= age[root]) return false;
    double dP = originPrior.ln(newO) - lnRoot; if (!isfinite(dP)) return false;
    double dF = P.lnD(newO) - P.lnD(oldO);
    if (log(rng.u()) < dF + dP) { origin = newO; lnFBD += dF; lnRoot += dP; return true; }
    return false;
}

bool Chain::moveRate() {
    if (K == 0) return false;
    int b = rng.below(K);
    double c = exp(tRate.val * (rng.u() - 0.5)), oldR = rates[b], newR = oldR * c;
    if (newR < 1e-100 || newR > 1e100) return false;
    double dP = lnRateP1(newR) - lnRateP1(oldR), dD = 0;
    rates[b] = newR;
    int term = bTerm[b]; double nv = 0;
    if (useData && term >= 0) { nv = termLL(term); dD = nv - tll[term]; }
    if (log(rng.u()) < dP + dD + log(c)) {
        lnRatePrior += dP; lnData += dD; sumRates += newR - oldR; if (useData && term >= 0) tll[term] = nv; return true;
    }
    rates[b] = oldR; return false;
}

bool Chain::moveMean() {
    if (K == 0) return false;
    double c = exp(tMean.val * (rng.u() - 0.5)), oldM = m, newM = m * c;
    if (newM < 1e-100 || newM > 1e100) return false;
    double nP = rateLnPriorGiven(rates, newM, sigma);
    double dP = nP - lnRatePrior;
    if (log(rng.u()) < dP + log(c)) { m = newM; uclnMu = log(m) - 0.5 * sigma * sigma; lnRatePrior = nP; return true; }
    return false;
}

bool Chain::moveSigma() {
    if (K == 0 || !ucln) return false;
    double c = exp(tSigma.val * (rng.u() - 0.5)), oldS = sigma, newS = sigma * c;
    double nP = rateLnPriorGiven(rates, m, newS);
    double dP = nP - lnRatePrior;
    if (log(rng.u()) < dP + log(c)) { sigma = newS; uclnMu = log(m) - 0.5 * sigma * sigma; lnRatePrior = nP; return true; }
    return false;
}

bool Chain::moveJoint() {
    if (K == 0) return false;
    double c = exp(tJoint.val * (rng.u() - 0.5)), oldM = m, newM = m * c;
    if (newM < 1e-100 || newM > 1e100) return false;
    vector<double> old = rates;
    for (auto &r : rates) r *= c;
    double nP = rateLnPriorGiven(rates, newM, sigma);
    double dP = nP - lnRatePrior, dD = 0; vector<double> nv;
    if (useData) { nv.resize(J); for (int j = 0; j < J; j++) { nv[j] = termLL(j); dD += nv[j] - tll[j]; } }
    double hast = (K + 1) * log(c);
    if (log(rng.u()) < dP + dD + hast) {
        m = newM; uclnMu = log(m) - 0.5 * sigma * sigma; lnRatePrior = nP; lnData += dD; sumRates *= c; if (useData) tll = nv; return true;
    }
    rates = old; return false;
}

bool Chain::moveVector() {
    if (K == 0) return false;
    double c = exp(tVec.val * (rng.u() - 0.5));
    vector<double> old = rates;
    for (auto &r : rates) r *= c;
    double nP = rateLnPriorGiven(rates, m, sigma);
    double dP = nP - lnRatePrior, dD = 0; vector<double> nv;
    if (useData) { nv.resize(J); for (int j = 0; j < J; j++) { nv[j] = termLL(j); dD += nv[j] - tll[j]; } }
    double hast = K * log(c);
    if (log(rng.u()) < dP + dD + hast) { lnRatePrior = nP; lnData += dD; sumRates *= c; if (useData) tll = nv; return true; }
    rates = old; return false;
}

bool Chain::moveFBDScale(int which) {
    Tune &t = which == 0 ? tLam : which == 1 ? tMu : tPsi;
    double &par_ = which == 0 ? lam : which == 1 ? mu : psi;
    double c = exp(t.val * (rng.u() - 0.5)), oldV = par_, newV = oldV * c;
    double dH = -10.0 * (newV - oldV);
    par_ = newV;
    FBD p; p.set(lam, mu, psi, rho);
    double nF = lnFBDFull(p);
    if (log(rng.u()) < (nF - lnFBD) + dH + log(c)) { P = p; lnFBD = nF; lnHyper += dH; return true; }
    par_ = oldV; return false;
}

bool Chain::moveUpDown() {
    double c = exp(tUD.val * (rng.u() - 0.5)), oLam = lam, oMu = mu, nLam = lam * c, nMu = mu * c;
    double dH = -10.0 * (nLam - oLam) - 10.0 * (nMu - oMu);
    lam = nLam; mu = nMu;
    FBD p; p.set(lam, mu, psi, rho);
    double nF = lnFBDFull(p);
    if (log(rng.u()) < (nF - lnFBD) + dH + 2.0 * log(c)) { P = p; lnFBD = nF; lnHyper += dH; return true; }
    lam = oLam; mu = oMu; return false;
}

void Chain::setupMoves() {
    weights.assign(N_MOVES, 0.0);
    double ni = N - 1, nf = mF, ne = (double)extNode.size();
    weights[M_NODE] = ni;
    if (hasClock) weights[M_COUPLED] = ne;
    weights[M_FOSAGE] = nf;
    weights[M_REGRAFT_L] = nf;
    weights[M_REGRAFT_G] = 2.0 * nf;
    weights[M_COLLEXP] = max(1.0, nf / 2.0);
    if (!useOrigin) { weights[M_ROOT] = 5 + ne / 50; if (hasClock) weights[M_ROOTC] = 5 + ne / 50; } else weights[M_ORIGIN] = 5;
    if (opt.fixFbd) { weights[M_ROOT] = weights[M_ROOTC] = weights[M_ORIGIN] = 0; }
    if (hasClock) { weights[M_RATE] = K; weights[M_RATESCALE_JOINT] = 3; weights[M_RATESCALE_ALL] = 2; weights[M_MEAN] = 5; if (ucln) weights[M_SIGMA] = 5; }
    if (!opt.fixFbd) { weights[M_LAMBDA] = 3; weights[M_MU] = 3; weights[M_PSI] = 3; weights[M_UPDOWN] = 2; }
    for (int i = 0; i < N_MOVES; i++) if (opt.noMoves.find(MOVE_NAMES[i]) != string::npos) weights[i] = 0.0;
    cum.resize(N_MOVES); totalW = 0; for (int i = 0; i < N_MOVES; i++) { totalW += weights[i]; cum[i] = totalW; }
}

void Chain::validate() {
    // structural + density consistency check (slow): tree topology valid, labels/regions/constraint nodes correct, caches correct
    // ages and SA flags
    int nsa = 0;
    if (par[root] != -1) die("root has a parent");
    for (int v = 0; v < nNodes; v++) {
        if (v != root) { int p = par[v]; if (p < N || (c0[p] != v && c1[p] != v)) die("broken parent link"); if (age[v] > age[p] + 1e-12) { char b[200]; snprintf(b, 200, "age order violated at %d: %.10g > %.10g", v, age[v], age[p]); die(b); } }
        if (v >= N) { if (c0[v] == c1[v]) die("twin children"); if (par[c0[v]] != v || par[c1[v]] != v) die("child link"); }
    }
    for (int f = n; f < N; f++) {
        if (age[f] < winMin[f] - 1e-12 || age[f] > winMax[f] + 1e-12) die("fossil outside window " + names[f]);
        if (sa[f]) { nsa++; if (age[par[f]] != age[f]) die("SA age mismatch"); if (isSAParent(sib(f))) {} if (par[f] == root && !useOrigin) die("root SA"); }
        else if (age[par[f]] <= age[f]) { /* allowed only if numerically equal */ }
    }
    for (int v = N; v < nNodes; v++) if (c0[v] < N && c1[v] < N && sa[c0[v]] && sa[c1[v]]) die("double SA");
    // constraint sets
    vector<vector<int>> tips(nNodes);
    function<void(int)> rec;
    vector<int> order; order.reserve(nNodes);
    {
        vector<int> st{root}; while (!st.empty()) { int v = st.back(); st.pop_back(); order.push_back(v); if (v >= N) { st.push_back(c0[v]); st.push_back(c1[v]); } }
    }
    if ((int)order.size() != nNodes) die("tree not connected: " + to_string(order.size()));
    for (auto it = order.rbegin(); it != order.rend(); ++it) { int v = *it; if (v < N) tips[v] = {v}; else { tips[v] = tips[c0[v]]; tips[v].insert(tips[v].end(), tips[c1[v]].begin(), tips[c1[v]].end()); sort(tips[v].begin(), tips[v].end()); } }
    vector<vector<char>> mem(nCons + 1, vector<char>(N, 0));
    for (int k = 0; k <= nCons; k++) for (int t : consMembers[k]) mem[k][t] = 1;
    vector<int> newLabel(nNodes, -1); vector<int> newX(nCons + 1, -1);
    for (int k = 1; k <= nCons; k++) {
        vector<int> s = consMembers[k]; sort(s.begin(), s.end());
        for (int v = 0; v < nNodes; v++) if (tips[v] == s) { newX[k] = v; break; }
        if (newX[k] < 0) die("constraint violated: " + consName[k]);
        if (newX[k] != x[k]) die("cached constraint node wrong " + consName[k]);
    }
    for (int v = 0; v < nNodes; v++) {
        if (v == root && !useOrigin) { if (label[v] != -1) die("root label"); continue; }
        int L = 0; size_t bs = N + 1;
        for (int k = 1; k <= nCons; k++) if (consMembers[k].size() < bs && consMembers[k].size() > tips[v].size()) {
            bool sub = true; for (int q : tips[v]) if (!mem[k][q]) { sub = false; break; }
            if (sub) { L = k; bs = consMembers[k].size(); }
        }
        for (int k = 1; k <= nCons; k++) if (x[k] == v) L = parentCons[k];
        if (L != label[v]) { char b[200]; snprintf(b, 200, "label wrong at node %d: have %d want %d", v, label[v], L); die(b); }
        if (region[L].empty() || pos[v] < 0 || region[L][pos[v]] != v) die("region array inconsistent");
    }
    // extant topology: each extant internal node keeps exactly its clade of extant tips with extant on both sides
    vector<int> extCnt(nNodes, 0);
    for (auto it = order.rbegin(); it != order.rend(); ++it) { int v = *it; if (v < N) extCnt[v] = v < n ? 1 : 0; else extCnt[v] = extCnt[c0[v]] + extCnt[c1[v]]; }
    for (int v : extNode) if (extCnt[c0[v]] == 0 || extCnt[c1[v]] == 0) die("extant node lost an extant child");
    int ne = 0; for (int v = N; v < nNodes; v++) if (extCnt[c0[v]] > 0 && extCnt[c1[v]] > 0) { ne++; if (!isExtInt[v]) die("non-extant node became extant"); }
    if (ne != (int)extNode.size()) die("extant internal node count changed");
    if (nsa == 0) {}
    // caches
    double a = lnFBDFull(P), r = useOrigin ? originPrior.ln(origin) : rootPrior.ln(age[root]);
    if (fabs(a - lnFBD) > 1e-6) { char b[300]; snprintf(b, 300, "lnFBD cache %.12g vs full %.12g", lnFBD, a); die(b); }
    if (fabs(r - lnRoot) > 1e-6) die("lnRoot cache");
    if (fabs(lnHyper - lnHyperFull()) > 1e-6) die("lnHyper cache");
    if (hasClock && fabs(rateLnPrior() - lnRatePrior) > 1e-6) die("rate prior cache");
    double d = 0; for (int j = 0; j < J; j++) { double t = useData ? termLL(j) : 0.0; d += t; if (fabs(t - tll[j]) > 1e-6) die("term cache"); } if (fabs(d - lnData) > 1e-6) die("data cache");
}

// ---- output -----------------------------------------------------------------------------------
string Chain::newick() const {
    string s; s.reserve(nNodes * 28);
    char buf[64];
    function<void(int)> rec = [&](int v) {
        if (v < N) { s += names[v]; } else { s += '('; rec(c0[v]); s += ','; rec(c1[v]); s += ')'; }
        double bl = (v == root) ? 0.0 : age[par[v]] - age[v];
        if (v != root) { if (v < N && sa[v]) s += ":0.0"; else { snprintf(buf, 64, ":%.12g", bl); s += buf; } }
    };
    rec(root); s += ';';
    return s;
}

void Chain::writeHeader(FILE *lf, FILE *tf) {
    fprintf(lf, "Iteration\tPosterior\tLikelihood\tPrior\tcond_age\troot_age_tree\tnum_samp_anc\tspeciation_rate\textinction_rate\tpsi\tdiversification\tturnover");
    if (useOrigin) fprintf(lf, "\torigin_time");
    for (auto &nm : extName) fprintf(lf, "\tage_%s", nm.c_str());
    if (hasClock) { fprintf(lf, "\tbranch_rates_mean\tmean_rt\tmean_rate"); if (ucln) fprintf(lf, "\tucln_sigma"); if (opt.logRates) for (int b = 0; b < K; b++) fprintf(lf, "\tbranch_rates[%d]", b + 1); }
    if (opt.logFossils) for (int f = n; f < N; f++) fprintf(lf, "\tage_%s\tsa_%s", names[f].c_str(), names[f].c_str());
    fprintf(lf, "\n");
    fprintf(tf, "Iteration\tPosterior\tLikelihood\tPrior\tfbd_tree\n");
}

void Chain::writeLog(FILE *lf) {
    double like = lnWin + lnData, prior = lnFBD + lnRoot + lnHyper + lnRatePrior;
    int nsa = 0; for (int f = n; f < N; f++) nsa += sa[f];
    double cond = useOrigin ? origin : age[root];
    double meanR = K > 0 ? accumulate(rates.begin(), rates.end(), 0.0) / K : 0.0;
    fprintf(lf, "%ld\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%d\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g", iter, like + prior, like, prior, cond, age[root], nsa, lam, mu, psi, lam - mu, mu / lam);
    if (useOrigin) fprintf(lf, "\t%.10g", origin);
    for (int v : extNode) fprintf(lf, "\t%.10g", age[v]);
    if (hasClock) { fprintf(lf, "\t%.10g\t%.10g\t%.10g", m, meanR, meanR); if (ucln) fprintf(lf, "\t%.10g", sigma); if (opt.logRates) for (int b = 0; b < K; b++) fprintf(lf, "\t%.10g", rates[b]); }
    if (opt.logFossils) for (int f = n; f < N; f++) fprintf(lf, "\t%.10g\t%d", age[f], (int)sa[f]);
    fprintf(lf, "\n");
}

void Chain::writeTree(FILE *tf) {
    double like = lnWin + lnData, prior = lnFBD + lnRoot + lnHyper + lnRatePrior;
    string nw = newick();
    fprintf(tf, "%ld\t%.10g\t%.10g\t%.10g\t%s\n", iter, like + prior, like, prior, nw.c_str());
}

void Chain::dumpState(FILE *f) {
    double cond = useOrigin ? origin : age[root];
    fprintf(f, "%ld\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t", iter, lam, mu, psi, cond, m);
    for (int b = 0; b < K; b++) fprintf(f, b ? ",%.17g" : "%.17g", rates[b]);
    double fb = lnFBDFull(P), r = useOrigin ? originPrior.ln(origin) : rootPrior.ln(age[root]);
    double d = 0; for (int j = 0; j < J; j++) d += termLL(j); // data terms reported even under --no-data, for the density check
    double rp = hasClock ? rateLnPrior() : 0.0;
    int nsa = 0; for (int q = n; q < N; q++) nsa += sa[q];
    fprintf(f, "\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%d\t", fb, r, lnHyperFull(), rp, d, lnWin, nsa);
    // full precision newick
    string s; char buf[64];
    function<void(int)> rec = [&](int v) {
        if (v < N) s += names[v]; else { s += '('; rec(c0[v]); s += ','; rec(c1[v]); s += ')'; }
        if (v != root) { double bl = age[par[v]] - age[v]; if (v < N && sa[v]) s += ":0.0"; else { snprintf(buf, 64, ":%.17g", bl); s += buf; } }
    };
    rec(root); s += ';';
    fprintf(f, "%s\t%.17g\n", s.c_str(), sigma);
}

void Chain::saveCheckpoint(const string &path) const {
    string tmp = path + ".tmp";
    {
        ofstream o(tmp); o << setprecision(17);
        o << iter << ' ' << burnDone << '\n' << lam << ' ' << mu << ' ' << psi << ' ' << m << ' ' << sigma << ' ' << origin << ' ' << root << '\n';
        o << rng.e << '\n';
        for (int v = 0; v < nNodes; v++) o << par[v] << ' ' << c0[v] << ' ' << c1[v] << ' ' << age[v] << ' ' << int(sa[v]) << ' ' << label[v] << '\n';
        for (int k = 0; k <= nCons; k++) o << x[k] << ' ';
        o << '\n';
        for (int k = 0; k <= nCons; k++) { o << region[k].size(); for (int v : region[k]) o << ' ' << v; o << '\n'; }
        for (double r : rates) o << r << ' '; o << '\n';
        const Tune *ts[] = {&tLam, &tMu, &tPsi, &tUD, &tRootS, &tRootCo, &tOrigin, &tRate, &tMean, &tJoint, &tVec, &tSigma};
        for (auto *t : ts) o << t->val << ' ' << t->att << ' ' << t->acc << ' ' << t->totalAtt << ' ' << t->totalAcc << ' ' << t->cycles << '\n';
        for (int i = 0; i < N_MOVES; i++) o << att[i] << ' ' << acc[i] << ' '; o << '\n';
        o << "END\n";
    }
    rename(tmp.c_str(), path.c_str());
}

void Chain::loadCheckpoint(const string &path) {
    ifstream in(path); if (!in) die("cannot read checkpoint " + path);
    in >> iter >> burnDone >> lam >> mu >> psi >> m >> sigma >> origin >> root >> rng.e;
    uclnMu = log(m) - 0.5 * sigma * sigma;
    for (int v = 0; v < nNodes; v++) { int s, l; in >> par[v] >> c0[v] >> c1[v] >> age[v] >> s >> l; sa[v] = (char)s; label[v] = -2 - l; }
    for (int v = 0; v < nNodes; v++) { label[v] = -2 - label[v]; pos[v] = -1; }
    for (int k = 0; k <= nCons; k++) in >> x[k];
    for (auto &r : region) r.clear();
    for (int k = 0; k <= nCons; k++) { size_t cnt; in >> cnt; region[k].resize(cnt); for (size_t i = 0; i < cnt; i++) { in >> region[k][i]; pos[region[k][i]] = (int)i; } }
    for (double &r : rates) in >> r;
    Tune *ts[] = {&tLam, &tMu, &tPsi, &tUD, &tRootS, &tRootCo, &tOrigin, &tRate, &tMean, &tJoint, &tVec, &tSigma};
    for (auto *t : ts) in >> t->val >> t->att >> t->acc >> t->totalAtt >> t->totalAcc >> t->cycles;
    for (int i = 0; i < N_MOVES; i++) in >> att[i] >> acc[i];
    string e; in >> e; if (e != "END") die("truncated checkpoint");
    recomputeAll();
}

void Chain::summary() const {
    printf("MOVE SUMMARY (proposals / accepted / rate / tuning):\n");
    const Tune *tmap[N_MOVES] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &tRootS, &tRootCo, &tOrigin, &tRate, &tJoint, &tVec, &tMean, &tSigma, &tLam, &tMu, &tPsi, &tUD};
    for (int i = 0; i < N_MOVES; i++) {
        if (weights[i] <= 0) continue;
        printf("  %-24s %12ld %12ld %8.4f", MOVE_NAMES[i], att[i], acc[i], att[i] ? (double)acc[i] / att[i] : 0.0);
        if (tmap[i]) printf("  tune=%.4g", tmap[i]->val);
        printf("\n");
    }
}

static string nowStr() { auto t = chrono::system_clock::to_time_t(chrono::system_clock::now()); string s = ctime(&t); s.pop_back(); return s; }

void Chain::run() {
    string logp = opt.prefix + ".log", tp = opt.prefix + ".trees", ckp = opt.prefix + ".ckp";
    FILE *lf, *tf, *df = nullptr;
    bool res = opt.resume;
    if (res) {
        loadCheckpoint(ckp);
        // truncate log and trees to entries <= iter
        for (const string &p : {logp, tp}) {
            ifstream in(p); string line; vector<string> keep; bool hdr = true;
            while (getline(in, line)) { if (hdr) { keep.push_back(line); hdr = false; continue; } if (line.empty()) continue; long it = atol(line.c_str()); if (it <= iter) keep.push_back(line); }
            in.close(); ofstream o(p); for (auto &l : keep) o << l << "\n";
        }
        lf = fopen(logp.c_str(), "a"); tf = fopen(tp.c_str(), "a");
        printf("RESUME from iteration %ld (%s)\n", iter, nowStr().c_str());
    } else {
        lf = fopen(logp.c_str(), "w"); tf = fopen(tp.c_str(), "w");
        writeHeader(lf, tf);
    }
    if (!opt.dumpFile.empty()) { df = fopen(opt.dumpFile.c_str(), res ? "a" : "w"); }
    if (!lf || !tf) die("cannot open output");
    printf("ffbd: taxa %d (extant %d, fossils %d), mode %s, cond %s, data %d, rho %g, moves total weight %g\n", N, n, mF, mode.c_str(), cond.c_str(), (int)useData, rho, totalW);
    printf("START lnFBD %.8f lnRoot %.8f lnHyper %.8f lnRatePrior %.8f lnData %.8f lnWin %.8f total %.8f\n", lnFBD, lnRoot, lnHyper, lnRatePrior, lnData, lnWin, lnFBD + lnRoot + lnHyper + lnRatePrior + lnData + lnWin);
    if (opt.check) validate();
    auto t0 = chrono::steady_clock::now();
    long totalTarget = opt.gens;
    // burn-in with tuning (only on a fresh chain)
    if (!res && opt.burnin > 0) {
        tuning = true;
        for (long b = 0; b < opt.burnin; b++) {
            double r = rng.u() * totalW; int mv = int(lower_bound(cum.begin(), cum.end(), r) - cum.begin()); if (mv >= N_MOVES) mv = N_MOVES - 1;
            bool ok = false;
            switch (mv) {
                case M_NODE: ok = moveNode(); break; case M_COUPLED: ok = moveCoupled(); break; case M_FOSAGE: ok = moveFossilAge(); break;
                case M_REGRAFT_L: ok = moveRegraft(true); break; case M_REGRAFT_G: ok = moveRegraft(false); break; case M_COLLEXP: ok = moveCollapseExpand(); break;
                case M_ROOT: ok = moveRoot(false); break; case M_ROOTC: ok = moveRoot(true); break; case M_ORIGIN: ok = moveOrigin(); break;
                case M_RATE: ok = moveRate(); break; case M_RATESCALE_JOINT: ok = moveJoint(); break; case M_RATESCALE_ALL: ok = moveVector(); break; case M_MEAN: ok = moveMean(); break; case M_SIGMA: ok = moveSigma(); break;
                case M_LAMBDA: ok = moveFBDScale(0); break; case M_MU: ok = moveFBDScale(1); break; case M_PSI: ok = moveFBDScale(2); break; case M_UPDOWN: ok = moveUpDown(); break;
            }
            att[mv]++; if (ok) acc[mv]++;
            Tune *tm[N_MOVES] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &tRootS, &tRootCo, &tOrigin, &tRate, &tJoint, &tVec, &tMean, &tSigma, &tLam, &tMu, &tPsi, &tUD};
            if (tm[mv]) { bool slide = (mv == M_ROOT || mv == M_ROOTC || mv == M_ORIGIN); tm[mv]->hit(ok); tm[mv]->adapt(slide ? 0.44 : 0.30, 1e-6, slide ? 200.0 : 2.0); }
            if (opt.checkPerStep) validate();
        }
        tuning = false;
        fill(att.begin(), att.end(), 0L); fill(acc.begin(), acc.end(), 0L);
        burnDone = opt.burnin;
        recomputeAll();
        printf("BURNIN done (%ld proposals), lnPost %.6f\n", opt.burnin, lnFBD + lnRoot + lnHyper + lnRatePrior + lnData + lnWin);
    }
    if (!res) { writeLog(lf); writeTree(tf); }
    long start = iter;
    for (long it = start + 1; it <= totalTarget; it++) {
        iter = it;
        double r = rng.u() * totalW; int mv = int(lower_bound(cum.begin(), cum.end(), r) - cum.begin()); if (mv >= N_MOVES) mv = N_MOVES - 1;
        bool ok = false;
        switch (mv) {
            case M_NODE: ok = moveNode(); break; case M_COUPLED: ok = moveCoupled(); break; case M_FOSAGE: ok = moveFossilAge(); break;
            case M_REGRAFT_L: ok = moveRegraft(true); break; case M_REGRAFT_G: ok = moveRegraft(false); break; case M_COLLEXP: ok = moveCollapseExpand(); break;
            case M_ROOT: ok = moveRoot(false); break; case M_ROOTC: ok = moveRoot(true); break; case M_ORIGIN: ok = moveOrigin(); break;
            case M_RATE: ok = moveRate(); break; case M_RATESCALE_JOINT: ok = moveJoint(); break; case M_RATESCALE_ALL: ok = moveVector(); break; case M_MEAN: ok = moveMean(); break; case M_SIGMA: ok = moveSigma(); break;
            case M_LAMBDA: ok = moveFBDScale(0); break; case M_MU: ok = moveFBDScale(1); break; case M_PSI: ok = moveFBDScale(2); break; case M_UPDOWN: ok = moveUpDown(); break;
        }
        att[mv]++; if (ok) acc[mv]++;
        if (opt.checkPerStep || (opt.check && it % 1000 == 0)) validate();
        if (opt.printEvery > 0 && it % opt.printEvery == 0) writeLog(lf);
        if (opt.treeEvery > 0 && it % opt.treeEvery == 0) { writeTree(tf); fflush(tf); fflush(lf); }
        if (df && opt.dumpEvery > 0 && it % opt.dumpEvery == 0) { dumpState(df); fflush(df); }
        if (opt.ckpEvery > 0 && it % opt.ckpEvery == 0) { fflush(lf); fflush(tf); saveCheckpoint(ckp); }
        if (it % 5000000 == 0) { // numerical refresh: bound floating-point drift of cached sums
            double before = lnFBD + lnData + lnRatePrior + lnRoot + lnHyper;
            recomputeAll();
            double after = lnFBD + lnData + lnRatePrior + lnRoot + lnHyper;
            if (fabs(after - before) > 1e-3) fprintf(stderr, "WARNING drift %.3g at %ld\n", after - before, it);
        }
    }
    auto t1 = chrono::steady_clock::now();
    double secs = chrono::duration<double>(t1 - t0).count();
    fclose(lf); fclose(tf); if (df) fclose(df);
    saveCheckpoint(ckp);
    printf("RUN done: %ld generations, s = %.3f\n", totalTarget - start, secs);
    summary();
    printf("FINAL lnPost %.6f\n", lnFBD + lnRoot + lnHyper + lnRatePrior + lnData + lnWin);
}

static void usage() {
    cerr << "usage: ffbd run INPUT --prefix P --seed S --gens G [--print N] [--tree-print N] [--checkpoint N] [--burnin N] [--mode fast|control] [--no-data] [--resume] [--check|--check-every-step] [--rb-sa-move] [--log-rates] [--dump FILE --dump-every N]\n";
    exit(2);
}

int main(int argc, char **argv) {
    if (argc < 3) usage();
    string cmd = argv[1];
    if (cmd != "run") usage();
    Options o; o.input = argv[2];
    for (int i = 3; i < argc; i++) {
        string a = argv[i];
        auto val = [&]() { if (i + 1 >= argc) usage(); return string(argv[++i]); };
        if (a == "--prefix") o.prefix = val(); else if (a == "--seed") o.seed = stoull(val()); else if (a == "--gens") o.gens = stol(val());
        else if (a == "--print") o.printEvery = stol(val()); else if (a == "--tree-print") o.treeEvery = stol(val()); else if (a == "--checkpoint") o.ckpEvery = stol(val());
        else if (a == "--burnin") o.burnin = stol(val()); else if (a == "--mode") o.mode = val(); else if (a == "--no-data") o.noData = true;
        else if (a == "--resume") o.resume = true; else if (a == "--check") o.check = true; else if (a == "--check-every-step") { o.check = true; o.checkPerStep = true; }
        else if (a == "--rb-sa-move") o.rbSA = true; else if (a == "--log-rates") o.logRates = true; else if (a == "--log-fossils") o.logFossils = true; else if (a == "--fix-fbd") o.fixFbd = true;
        else if (a == "--no-moves") o.noMoves = val(); else if (a == "--dump") o.dumpFile = val(); else if (a == "--dump-every") o.dumpEvery = stol(val());
        else usage();
    }
    Chain c(o); c.read(o.input); c.init(); c.run();
    return 0;
}
