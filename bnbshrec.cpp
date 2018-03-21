/*
 * A simple multithreaded interval-based bnb solver
 */

/* 
 * File:   tutorialbnb.cpp
 * Author: mposypkin
 *
 * Created on December 13, 2017, 3:22 PM
 */

#include <iostream>
#include <limits>
#include <random>
#include <algorithm>
#include <vector>
#include <iterator>
#include <functional>
#include <thread>
#include <chrono>
#include <shared_mutex>
#include <testfuncs/benchmarks.hpp>

using BM = Benchmark<double>;
using Box = std::vector<Interval<double>>;

static int procs = 64;

static int mtStepsLimit = 1000;

static int maxStepsTotal = 1000000;

struct Record {

    double getValue() const {
        //std::shared_lock<std::shared_mutex> lock(mMut);
        //std::unique_lock<std::shared_mutex> lock(mMut);
        std::lock_guard<std::mutex> lock(mMut);
        return mValue;
    }

    void update(double nval, const std::vector<double>& nrecord) {
        //std::unique_lock<std::shared_mutex> lock(mMut);
        std::lock_guard<std::mutex> lock(mMut);
        mValue = nval;
        mVector = nrecord;
    }

    //mutable std::shared_mutex mMut;
    mutable std::mutex mMut;
    double mValue;
    std::vector<double> mVector;
};

struct State {

    void merge(const State& s) {
        mSteps += s.mSteps;
        mPool.insert(mPool.end(), s.mPool.begin(), s.mPool.end());
    }

    void split(State& s1, State& s2) {
        const int remMaxSteps = mMaxSteps - mSteps;
        s1.mMaxSteps = remMaxSteps / 2;
        s2.mMaxSteps = remMaxSteps - s1.mMaxSteps;
        s1.mProcs = mProcs / 2;
        s2.mProcs = mProcs - s1.mProcs;

        while (true) {
            if (mPool.empty())
                break;
            s1.mPool.push_back(mPool.back());
            mPool.pop_back();
            if (mPool.empty())
                break;
            s2.mPool.push_back(mPool.back());
            mPool.pop_back();
        }
    }


    std::vector<Box> mPool;

    int mMaxSteps;

    int mSteps = 0;

    int mProcs;
};

Record record;

std::ostream& operator<<(std::ostream & out, const State s) {
    out << "\"steps\" :" << s.mSteps << "\n";
    out << "\"max steps\" :" << s.mMaxSteps << "\n";
    return out;
}

double len(const Interval<double>& I) {
    return I.rb() - I.lb();
}

void split(const Box& ibox, std::vector<Box>& v) {
    auto result = std::max_element(ibox.begin(), ibox.end(),
            [](const Interval<double>& f, const Interval<double>& s) {
                return len(f) < len(s);
            });
    const int i = result - ibox.begin();
    const double maxlen = len(ibox[i]);
    Box b1(ibox);
    Interval<double> ilow(ibox[i].lb(), ibox[i].lb() + 0.5 * maxlen);
    b1[i] = ilow;
    Box b2(ibox);
    Interval<double> iupper(ibox[i].lb() + 0.5 * maxlen, ibox[i].rb());
    b2[i] = iupper;
    v.push_back(std::move(b1));
    v.push_back(std::move(b2));
}

void getCenter(const Box& ibox, std::vector<double>& c) {
    const int n = ibox.size();
    for (int i = 0; i < n; i++) {
        c[i] = 0.5 * (ibox[i].lb() + ibox[i].rb());
    }
}

void solveSerial(State& s, const BM& bm, double eps) {
    const int dim = bm.getDim();
    std::vector<double> c(dim);
    while (!s.mPool.empty()) {
        s.mSteps++;
        Box b = s.mPool.back();
        s.mPool.pop_back();
        getCenter(b, c);
        double v = bm.calcFunc(c);
        if (v < record.getValue()) {
            record.update(v, c);
        }
        auto lb = bm.calcInterval(b).lb();
        if (lb <= record.getValue() - eps) {
            split(b, s.mPool);
        }
        if (s.mSteps >= s.mMaxSteps)
            break;
    }
}

void solve(State& s, const BM& bm, double eps) {
    if ((s.mProcs == 1) || (s.mMaxSteps <= mtStepsLimit)) {
        solveSerial(s, bm, eps);
    } else {
        auto presolve = [&](State & s) {
            const int dim = bm.getDim();
            std::vector<double> c(dim);
            while ((!s.mPool.empty()) && (s.mPool.size() < 2)) {
                s.mSteps++;
                Box b = s.mPool.back();
                s.mPool.pop_back();
                getCenter(b, c);
                double v = bm.calcFunc(c);
                if (v < record.getValue()) {
                    record.update(v, c);
                }
                auto lb = bm.calcInterval(b).lb();
                if (lb <= record.getValue() - eps) {
                    split(b, s.mPool);
                }
                if (s.mSteps >= s.mMaxSteps)
                    break;
            }
        };

        while (true) {
            if (s.mPool.empty())
                break;
            //std::cout << s << "\n";
            presolve(s);
            const int remMaxSteps = s.mMaxSteps - s.mSteps;
            if (remMaxSteps == 0)
                break;
            if (remMaxSteps <= mtStepsLimit) {
                solveSerial(s, bm, eps);
                break;
            } else {
                State s1, s2;
                s.split(s1, s2);

#if 1            
                std::thread t1(solve, std::ref(s1), std::ref(bm), eps);
                std::thread t2(solve, std::ref(s2), std::ref(bm), eps);
                t1.join();
                t2.join();
#else
                solve(s1, bm, eps);
                solve(s2, bm, eps);
#endif
                s.merge(s1);
                s.merge(s2);
            }
        }
    }
}

double findMin(const BM& bm, double eps, int maxstep) {
    const int dim = bm.getDim();
    Box ibox;
    for (int i = 0; i < dim; i++) {
        ibox.emplace_back(bm.getBounds()[i].first, bm.getBounds()[i].second);
    }
    State s;
    std::vector<double> recvec(dim, 0);
    record.update(std::numeric_limits<double>::max(), recvec);
    s.mPool.push_back(ibox);
    s.mMaxSteps = maxstep;
    s.mProcs = procs;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
#if 0   
    solveSerial(s, bm, eps);
#else    
    solve(s, bm, eps);
#endif
    end = std::chrono::system_clock::now();
    int mseconds = (std::chrono::duration_cast<std::chrono::microseconds> (end - start)).count();
    std::cout << "Time: " << mseconds << " microsecond\n";
    std::cout << "Time per subproblem: " << (double) s.mSteps / (double) mseconds << " miscroseconds." << std::endl;
    if (s.mSteps >= maxstep) {
        std::cout << "Failed to converge in " << maxstep << " steps\n";
    } else {
        std::cout << "Converged in " << s.mSteps << " steps\n";
    }

    std::cout << "BnB found = " << record.getValue() << std::endl;
    std::cout << " at x [ ";
    std::copy(record.mVector.begin(), record.mVector.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "]\n";
    return record.getValue();
}

bool testBench(const BM& bm) {
    constexpr double eps = 0.1;
    std::cout << "*************Testing benchmark**********" << std::endl;
    std::cout << bm;
    double v = findMin(bm, eps, maxStepsTotal);
    double diff = v - bm.getGlobMinY();
    if (diff > eps) {
        std::cout << "BnB failed for " << bm.getDesc() << " benchmark " << std::endl;
    }
    std::cout << "the difference is " << v - bm.getGlobMinY() << std::endl;
    std::cout << "****************************************" << std::endl << std::endl;
}

main(int argc, char* argv[]) {
    if (argc >= 2) {
        procs = atoi(argv[1]);
    }
    if (argc >= 3) {
        mtStepsLimit = atoi(argv[2]);
    }
    if (argc >= 4) {
        maxStepsTotal = atoi(argv[3]);
    }
    std::cout << "Smart solver with np = " << procs << ", mtStepsLimit =  " << mtStepsLimit << ", maxStepsTotal = " << maxStepsTotal << std::endl;
    
#if 0
    PowellSingular2Benchmark<double> pb(8);
    testBench(pb);
#else    
    Benchmarks<double> tests;
    for (auto bm : tests) {
        testBench(*bm);
    }
#endif  
    
}
