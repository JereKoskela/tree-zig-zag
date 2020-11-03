#ifndef SIM
#define SIM

#include <algorithm>
#include <cassert>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_spmatrix.h>
#include <iostream>
#include <limits>
#include <sstream>
#include <unistd.h>
#include <vector>

struct Sim {

    Sim(std::string filename, const double mutation_v, const double hyb,
        const double mutation_sd)
        : n(0), nsites(), hybrid_rate(hyb), mutation_rate(),
          mutation_vel(mutation_v), hybrid_mutation_sd(mutation_sd),
          row_count(), left_child(), right_child(), t(), v(), l(), L(), m(),
          M(), L_tmp(), M_tmp() {
        gen = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(gen, time(NULL) * getpid());
        data = gsl_spmatrix_alloc(1, 1);
        std::string line, token;
        int s = 0;
        std::ifstream file;
        file.open(filename);
        while (getline(file, line)) {
            std::vector<int> row;
            std::stringstream iss;
            iss << line;
            while (getline(iss, token, ' ')) {
                row.push_back(atoi(token.c_str()));
            }
            for (unsigned int i = 0; i < row.size() - 1; i++) {
                if (row[i] == 1) {
                    gsl_spmatrix_set(data, s, i, 1);
                }
            }
            row_count.push_back(row.back());
            n += row.back();
            s++;
        }
        nsites = data->size2;
        left_child.resize(n - 1, -1);
        right_child.resize(n - 1, -1);
        v.resize(n - 1, 1);
        for (int i = 0; i < n - 1; i++) {
            v[i] = 1 / gsl_sf_choose(n - i, 2);
            if (gsl_rng_uniform(gen) < 0.5) {
                v[i] *= -1;
            }
        }
        mutation_vel /= nsites;
        hybrid_mutation_sd /= nsites;
        if (gsl_rng_uniform(gen) < 0.5) {
            mutation_vel *= -1;
        }
        t.resize(n - 1, 0);
        std::vector<double> tmp(2);
        std::vector<std::vector<double>> tmp_2(2 * n - 1, tmp);
        L_tmp.resize(nsites, tmp_2);
        M_tmp.resize(nsites, tmp_2);
        l.resize(nsites, tmp_2);
        L.resize(nsites, tmp_2);
        m.resize(nsites, tmp_2);
        M.resize(nsites, tmp_2);
        for (int i = 0; i < nsites; i++) {
            for (int j = 0; j < n; j++) {
                if (gsl_spmatrix_get(data, leaf_to_row(j), i) == 0) {
                    L[i][j][0] = 1;
                    L[i][j][1] = 0;
                    l[i][j][0] = 1;
                    l[i][j][1] = 0;
                    L_tmp[i][j][0] = 1;
                    L_tmp[i][j][1] = 0;
                } else {
                    L[i][j][0] = 0;
                    L[i][j][1] = 1;
                    l[i][j][0] = 0;
                    l[i][j][1] = 1;
                    L_tmp[i][j][0] = 0;
                    L_tmp[i][j][1] = 1;
                }
            }
            M_tmp[i][2 * n - 2][0] = 0.5;
            M_tmp[i][2 * n - 2][1] = 0.5;
            m[i][2 * n - 2][0] = 0.5;
            m[i][2 * n - 2][1] = 0.5;
            M[i][2 * n - 2][0] = 0.5;
            M[i][2 * n - 2][1] = 0.5;
        }
        mutation_rate = watterson_estimator() / nsites;
        sample_tree();
    }

    ~Sim() {
        gsl_rng_free(gen);
        gsl_spmatrix_free(data);
    }

    double watterson_estimator() const {
        double harmonic = 1;
        for (int i = 2; i < n; i++) {
            harmonic += 1.0 / i;
        }
        int seg = 0;
        int anc;
        for (int i = 0; i < nsites; i++) {
            anc = gsl_spmatrix_get(data, 0, i);
            for (unsigned int j = 1; j < row_count.size(); j++) {
                if (gsl_spmatrix_get(data, j, i) != anc) {
                    seg++;
                    break;
                }
            }
        }
        return seg / harmonic;
    }

    void sample_tree() {
        std::vector<int> active(n);
        for (int i = 0; i < n; i++) {
            active[i] = i;
        }
        int next_parent = n;
        for (int i = 0; i < (int)row_count.size(); i++) {
            for (int j = 0; j < row_count[i] - 1; j++) {
                left_child[next_parent - n] = active[i];
                right_child[next_parent - n] = active[i + 1];
                t[next_parent - n] = gsl_ran_exponential(
                    gen, 1 / gsl_sf_choose(active.size(), 2));
                active[i] = next_parent;
                active.erase(active.begin() + i + 1);
                next_parent++;
            }
        }
        for (int i = 0; i < (int)row_count.size() - 1; i++) {
            left_child[next_parent - n] = active[0];
            right_child[next_parent - n] = active[1];
            t[next_parent - n] =
                gsl_ran_exponential(gen, 1 / gsl_sf_choose(active.size(), 2));
            active[0] = next_parent;
            active.erase(active.begin() + 1);
            next_parent++;
        }
        return;
    }

    int leaf_to_row(const int leaf) const {
        int row = 0;
        int total = row_count[0];
        while (leaf + 1 > total) {
            row++;
            total += row_count[row];
        }
        return row;
    }

    double branch_length(const int i, const double add) const {
        int ind = std::max(i - n + 1, 0);
        double ret = t[ind] + v[ind] * add;
        while (left_child[ind] != i && right_child[ind] != i) {
            ind++;
            ret += t[ind] + v[ind] * add;
        }
        return ret;
    }

    double branch_velocity(const int i) const {
        int ind = std::max(i - n + 1, 0);
        double ret = v[ind];
        while (left_child[ind] != i && right_child[ind] != i) {
            ind++;
            ret += v[ind];
        }
        return ret;
    }

    double post_order(const int i, const int c, const int h,
                      const double s) const {
        double bl = branch_length(c, s);
        int g = (h + 1) % 2;
        double ret = (L_tmp[i][c][0] + L_tmp[i][c][1] +
                      (L_tmp[i][c][h] - L_tmp[i][c][g]) *
                          exp(-(mutation_rate + mutation_vel * s) * bl)) /
                     2;
        return ret;
    }

    double post_order_ub(const int i, const int c, const int h,
                         const double s) const {
        double ret = L[i][c][0] + L[i][c][1];
        double bl = branch_length(c, 0);
        double vl = branch_velocity(c);
        int g = (h + 1) % 2;
        if (L[i][c][h] > L[i][c][g]) {
            ret +=
                (L[i][c][h] - L[i][c][g]) *
                exp(-fmin(mutation_rate * bl,
                          (mutation_rate + mutation_vel * s) * (bl + vl * s)));
        } else {
            double crit = -bl / (2 * vl) - mutation_rate / (2 * mutation_vel);
            if (crit > 0 && crit < s) {
                ret += (L[i][c][h] - L[i][c][g]) *
                       exp(-(mutation_rate - mutation_vel * bl / vl) *
                           (bl - vl * mutation_rate / mutation_vel) / 4);
            } else {
                ret += (L[i][c][h] - L[i][c][g]) *
                       exp(-fmax(mutation_rate * bl,
                                 (mutation_rate + mutation_vel * s) *
                                     (bl + vl * s)));
            }
        }
        return fmin(ret / 2, 1);
    }

    double post_order_lb(const int i, const int c, const int h,
                         const double s) const {
        double ret = l[i][c][0] + l[i][c][1];
        double bl = branch_length(c, 0);
        double vl = branch_velocity(c);
        int g = (h + 1) % 2;
        if (l[i][c][h] > l[i][c][g]) {
            double crit = -bl / (2 * vl) - mutation_rate / (2 * mutation_vel);
            if (crit > 0 && crit < s) {
                ret += (l[i][c][h] - l[i][c][g]) *
                       exp(-(mutation_rate - mutation_vel * bl / vl) *
                           (bl - vl * mutation_rate / mutation_vel) / 4);
            } else {
                ret += (l[i][c][h] - l[i][c][g]) *
                       exp(-fmax(mutation_rate * bl,
                                 (mutation_rate + mutation_vel * s) *
                                     (bl + vl * s)));
            }
        } else {
            ret +=
                (l[i][c][h] - l[i][c][g]) *
                exp(-fmin(mutation_rate * bl,
                          (mutation_rate + mutation_vel * s) * (bl + vl * s)));
        }
        return fmax(ret / 2, 0);
    }

    double pre_order(const int i, const int p, const int c, const int d,
                     const int h, const double s) const {
        double bl = branch_length(c, s);
        int g = (h + 1) % 2;
        double c0 = M_tmp[i][p][h] * post_order(i, d, h, s);
        double c1 = M_tmp[i][p][g] * post_order(i, d, g, s);
        double ret =
            (c0 + c1 +
             (c0 - c1) * exp(-(mutation_rate + mutation_vel * s) * bl)) /
            2;
        return ret;
    }

    double pre_order_ub(const int i, const int p, const int c, const int d,
                        const int h, const double s) const {
        double bl = branch_length(c, 0);
        double vl = branch_velocity(c);
        int g = (h + 1) % 2;
        double c0 = M[i][p][h] * post_order_ub(i, d, h, s);
        double c1 = M[i][p][g] * post_order_ub(i, d, g, s);
        double ret = c0 + c1;
        if (c0 > c1) {
            ret +=
                (c0 - c1) *
                exp(-fmin(mutation_rate * bl,
                          (mutation_rate + mutation_vel * s) * (bl + vl * s)));
        } else {
            double crit = -bl / (2 * vl) - mutation_rate / (2 * mutation_vel);
            if (crit > 0 && crit < s) {
                ret += (c0 - c1) *
                       exp(-(mutation_rate - mutation_vel * bl / vl) *
                           (bl - vl * mutation_rate / mutation_vel) / 4);
            } else {
                ret +=
                    (c0 - c1) * exp(-fmax(mutation_rate * bl,
                                          (mutation_rate + mutation_vel * s) *
                                              (bl + vl * s)));
            }
        }
        return fmin(ret / 2, 1);
    }

    double pre_order_lb(const int i, const int p, const int c, const int d,
                        const int h, const double s) const {
        double bl = branch_length(c, 0);
        double vl = branch_velocity(c);
        int g = (h + 1) % 2;
        double c0 = m[i][p][h] * post_order_lb(i, d, h, s);
        double c1 = m[i][p][g] * post_order_lb(i, d, g, s);
        double ret = c0 + c1;
        if (c0 > c1) {
            double crit = -bl / (2 * vl) - mutation_rate / (2 * mutation_vel);
            if (crit > 0 && crit < s) {
                ret += (c0 - c1) *
                       exp(-(mutation_rate - mutation_vel * bl / vl) *
                           (bl - vl * mutation_rate / mutation_vel) / 4);
            } else {
                ret +=
                    (c0 - c1) * exp(-fmax(mutation_rate * bl,
                                          (mutation_rate + mutation_vel * s) *
                                              (bl + vl * s)));
            }
        } else {
            ret +=
                (c0 - c1) *
                exp(-fmin(mutation_rate * bl,
                          (mutation_rate + mutation_vel * s) * (bl + vl * s)));
        }
        return fmax(ret / 2, 0);
    }

    void fill_bounds(const double s) {
        int c, d;
        for (int i = 0; i < nsites; i++) {
            for (int j = n; j < 2 * n - 1; j++) {
                L[i][j][0] = post_order_ub(i, left_child[j - n], 0, s) *
                             post_order_ub(i, right_child[j - n], 0, s);
                L[i][j][1] = post_order_ub(i, left_child[j - n], 1, s) *
                             post_order_ub(i, right_child[j - n], 1, s);
                l[i][j][0] = post_order_lb(i, left_child[j - n], 0, s) *
                             post_order_lb(i, right_child[j - n], 0, s);
                l[i][j][1] = post_order_lb(i, left_child[j - n], 1, s) *
                             post_order_lb(i, right_child[j - n], 1, s);
            }
            for (int j = 2 * n - 2; j >= n; j--) {
                c = left_child[j - n];
                d = right_child[j - n];
                M[i][c][0] = pre_order_ub(i, j, c, d, 0, s);
                M[i][c][1] = pre_order_ub(i, j, c, d, 1, s);
                M[i][d][0] = pre_order_ub(i, j, d, c, 0, s);
                M[i][d][1] = pre_order_ub(i, j, d, c, 1, s);
                m[i][c][0] = pre_order_lb(i, j, c, d, 0, s);
                m[i][c][1] = pre_order_lb(i, j, c, d, 1, s);
                m[i][d][0] = pre_order_lb(i, j, d, c, 0, s);
                m[i][d][1] = pre_order_lb(i, j, d, c, 1, s);
            }
        }
        return;
    }

    double grad_log_L_ub(const int ind, const double s,
                         const std::vector<int> &branches) {
        double ret = 0;
        double temp = 0;
        if (ind >= 0) {
            int c;
            for (int i = 0; i < nsites; i++) {
                temp = 0;
                for (unsigned int j = 0; j < branches.size(); j++) {
                    c = branches[j];
                    if (v[ind] > 0) {
                        if (l[i][c][1] < L[i][c][0]) {
                            temp += M[i][c][0] * (l[i][c][1] - L[i][c][0]);
                        } else {
                            temp += m[i][c][0] * (l[i][c][1] - L[i][c][0]);
                        }
                        if (l[i][c][0] < L[i][c][1]) {
                            temp += M[i][c][1] * (l[i][c][0] - L[i][c][1]);
                        } else {
                            temp += m[i][c][1] * (l[i][c][0] - L[i][c][1]);
                        }
                    } else {
                        if (L[i][c][1] > l[i][c][0]) {
                            temp += M[i][c][0] * (L[i][c][1] - l[i][c][0]);
                        } else {
                            temp += m[i][c][0] * (L[i][c][1] - l[i][c][0]);
                        }
                        if (L[i][c][0] > l[i][c][1]) {
                            temp += M[i][c][1] * (L[i][c][0] - l[i][c][1]);
                        } else {
                            temp += m[i][c][1] * (L[i][c][0] - l[i][c][1]);
                        }
                    }
                }
                if (temp * v[ind] < 0) {
                    ret += (mutation_rate + fmax(0, mutation_vel * s)) * temp /
                           (l[i][2 * n - 2][0] + l[i][2 * n - 2][1]);
                } else {
                    ret += (mutation_rate + fmin(0, mutation_vel * s)) * temp /
                           (L[i][2 * n - 2][0] + L[i][2 * n - 2][1]);
                }
            }
        } else {
            double bl, vl, ret_i;
            for (int i = 0; i < nsites; i++) {
                ret_i = 0;
                for (int j = 0; j < 2 * n - 2; j++) {
                    temp = 0;
                    if (mutation_vel > 0) {
                        if (l[i][j][1] < L[i][j][0]) {
                            temp += M[i][j][0] * (l[i][j][1] - L[i][j][0]);
                        } else {
                            temp += m[i][j][0] * (l[i][j][1] - L[i][j][0]);
                        }
                        if (l[i][j][0] < L[i][j][1]) {
                            temp += M[i][j][1] * (l[i][j][0] - L[i][j][1]);
                        } else {
                            temp += m[i][j][1] * (l[i][j][0] - L[i][j][1]);
                        }
                    } else {
                        if (L[i][j][1] > l[i][j][0]) {
                            temp += M[i][j][0] * (L[i][j][1] - l[i][j][0]);
                        } else {
                            temp += m[i][j][0] * (L[i][j][1] - l[i][j][0]);
                        }
                        if (L[i][j][0] > l[i][j][1]) {
                            temp += M[i][j][1] * (L[i][j][0] - l[i][j][1]);
                        } else {
                            temp += m[i][j][1] * (L[i][j][0] - l[i][j][1]);
                        }
                    }
                    bl = branch_length(j, 0);
                    vl = branch_velocity(j);
                    if (temp * mutation_vel < 0) {
                        ret_i += (bl + fmax(0, vl * s)) * temp;
                    } else {
                        ret_i += (bl + fmin(0, vl * s)) * temp;
                    }
                }
                if (ret_i * mutation_vel < 0) {
                    ret += ret_i / (l[i][2 * n - 2][0] + l[i][2 * n - 2][1]);
                } else {
                    ret += ret_i / (L[i][2 * n - 2][0] + L[i][2 * n - 2][1]);
                }
            }
        }
        return ret;
    }

    double grad_log_L(const int ind, const double s,
                      const std::vector<int> &branches) {
        double ret = 0;
        int c, d;
        for (int i = 0; i < nsites; i++) {
            for (int j = n; j < 2 * n - 1; j++) {
                L_tmp[i][j][0] = post_order(i, left_child[j - n], 0, s) *
                                 post_order(i, right_child[j - n], 0, s);
                L_tmp[i][j][1] = post_order(i, left_child[j - n], 1, s) *
                                 post_order(i, right_child[j - n], 1, s);
            }
            for (int j = 2 * n - 2; j >= n; j--) {
                c = left_child[j - n];
                d = right_child[j - n];
                M_tmp[i][c][0] = pre_order(i, j, c, d, 0, s);
                M_tmp[i][c][1] = pre_order(i, j, c, d, 1, s);
                M_tmp[i][d][0] = pre_order(i, j, d, c, 0, s);
                M_tmp[i][d][1] = pre_order(i, j, d, c, 1, s);
            }
            if (ind >= 0) {
                for (unsigned int j = 0; j < branches.size(); j++) {
                    c = branches[j];
                    ret += (mutation_rate + mutation_vel * s) *
                           (M_tmp[i][c][0] - M_tmp[i][c][1]) *
                           (L_tmp[i][c][1] - L_tmp[i][c][0]) /
                           (L_tmp[i][2 * n - 2][0] + L_tmp[i][2 * n - 2][1]);
                }
            } else {
                for (int j = 0; j < 2 * n - 2; j++) {
                    ret += branch_length(j, s) *
                           (M_tmp[i][j][0] - M_tmp[i][j][1]) *
                           (L_tmp[i][j][1] - L_tmp[i][j][0]) /
                           (L_tmp[i][2 * n - 2][0] + L_tmp[i][2 * n - 2][1]);
                }
            }
        }
        return ret;
    }

    double sample_increment(const int i, const double t_max,
                            std::vector<int> &branches) {
        double ret = 0;
        double rate = 0;
        if (i == n - 1) {
            rate = -mutation_vel * grad_log_L_ub(-1, t_max, branches);
        } else {
            rate = v[i] * (gsl_sf_choose(n - i, 2) -
                           grad_log_L_ub(i, t_max, branches));
        }
        double alpha;
        do {
            if (rate > 0) {
                ret += gsl_ran_exponential(gen, 1 / rate);
            } else {
                ret = t_max;
            }
            if (ret >= t_max) {
                ret = t_max;
                alpha = 1;
            } else {
                if (i == n - 1) {
                    alpha =
                        -mutation_vel * grad_log_L(-1, ret, branches) / rate;
                } else {
                    alpha = v[i] *
                            (gsl_sf_choose(n - i, 2) -
                             grad_log_L(i, ret, branches)) /
                            rate;
                }
            }
        } while (gsl_rng_uniform(gen) > alpha);
        return ret;
    }

    double log_likelihood_types() {
        double ret = 0;
        for (int i = 0; i < nsites; i++) {
            for (int j = n; j < 2 * n - 1; j++) {
                L_tmp[i][j][0] = post_order(i, left_child[j - n], 0, 0) *
                                 post_order(i, right_child[j - n], 0, 0);
                L_tmp[i][j][1] = post_order(i, left_child[j - n], 1, 0) *
                                 post_order(i, right_child[j - n], 1, 0);
            }
            ret +=
                log(L_tmp[i][2 * n - 2][0] + L_tmp[i][2 * n - 2][1]) - log(2);
        }
        return ret;
    }

    double log_likelihood_times() const {
        double ret = 0;
        for (int i = 0; i < n - 1; i++) {
            ret -= gsl_sf_choose(n - i, 2) * t[i];
        }
        return ret;
    }

    int detach_reattach(const int detach, int attach_above,
                        const double new_time) {
        // Returns the index of the branch above which to detach-reattach the
        // detached branch in order to reverse the move.
        int parent = std::max(detach - n, 0);
        while (left_child[parent] != detach && right_child[parent] != detach) {
            parent++;
        }
        int sibling = left_child[parent];
        if (detach == sibling) {
            sibling = right_child[parent];
        }
        if (parent < n - 2) {
            t[parent + 1] += t[parent];
        }
        t.erase(t.begin() + parent);
        left_child.erase(left_child.begin() + parent);
        right_child.erase(right_child.begin() + parent);
        for (int i = parent; i < n - 2; i++) {
            if (left_child[i] == n + parent) {
                left_child[i] = sibling;
            } else if (left_child[i] > n + parent) {
                left_child[i]--;
            }
            if (right_child[i] == n + parent) {
                right_child[i] = sibling;
            } else if (right_child[i] > n + parent) {
                right_child[i]--;
            }
        }
        if (attach_above > n + parent) {
            attach_above--;
        }
        int ind = 0;
        double cumulative_time = t[0];
        while (ind < n - 2 && new_time > cumulative_time) {
            ind++;
            if (ind < n - 2) {
                cumulative_time += t[ind];
            }
        }
        if (ind < n - 2) {
            t.insert(t.begin() + ind, new_time - cumulative_time + t[ind]);
            t[ind + 1] = cumulative_time - new_time;
        } else {
            if (n == 2) {
                t.insert(t.end(), new_time);
            } else {
                t.insert(t.end(), new_time - cumulative_time);
            }
        }
        left_child.insert(left_child.begin() + ind, detach);
        right_child.insert(right_child.begin() + ind, attach_above);
        for (int i = ind + 1; i < n - 1; i++) {
            if (left_child[i] == attach_above) {
                left_child[i] = n + ind;
            } else if (left_child[i] >= n + ind) {
                left_child[i]++;
            }
            if (right_child[i] == attach_above) {
                right_child[i] = n + ind;
            } else if (right_child[i] >= n + ind) {
                right_child[i]++;
            }
        }
        if (attach_above >= n + ind) {
            attach_above++;
        }
        if (sibling >= n + ind) {
            sibling++;
        }
        return sibling;
    }

    int metropolis_spr() {
        int ret = 1;
        int parent = 0;
        if (n > 2) {
            parent = gsl_rng_uniform_int(gen, n - 2);
        }
        int detach, sibling;
        if (gsl_rng_uniform(gen) < 0.5) {
            sibling = left_child[parent];
            detach = right_child[parent];
        } else {
            sibling = right_child[parent];
            detach = left_child[parent];
        }
        int new_parent, attach_above;
        do {
            new_parent = -1;
            attach_above = 2 * n - 2;
            if (gsl_rng_uniform(gen) > 1.0 / (2 * n - 1)) {
                if (n == 2) {
                    new_parent = 0;
                } else {
                    new_parent = gsl_rng_uniform_int(gen, n - 2);
                }
                if (gsl_rng_uniform(gen) < 0.5) {
                    attach_above = left_child[new_parent];
                } else {
                    attach_above = right_child[new_parent];
                }
                while (attach_above == detach) {
                    if (n == 2) {
                        new_parent = 0;
                    } else {
                        new_parent = gsl_rng_uniform_int(gen, n - 2);
                    }
                    if (gsl_rng_uniform(gen) < 0.5) {
                        attach_above = left_child[new_parent];
                    } else {
                        attach_above = right_child[new_parent];
                    }
                }
            }
        } while (attach_above == n + parent);
        if (new_parent != -1 && detach - n >= new_parent) {
            ret = 0;
        } else {
            double log_alpha = -log_likelihood_times() - log_likelihood_types();
            double old_parent_time = 0;
            for (int i = 0; i <= parent; i++) {
                old_parent_time += t[i];
            }
            double new_parent_time = 0;
            for (int i = 0; i <= std::max(attach_above, detach) - n; i++) {
                new_parent_time += t[i];
            }
            double increment = 0;
            if (new_parent == -1) {
                increment = gsl_ran_exponential(gen, 1);
                log_alpha += increment;
            } else {
                double ub = new_parent_time;
                for (int i =
                         std::max(std::max(attach_above, detach) - n, -1) + 1;
                     i <= new_parent; i++) {
                    ub += t[i];
                }
                increment =
                    gsl_ran_flat(gen, new_parent_time, ub) - new_parent_time;
                log_alpha += log(ub - new_parent_time);
            }
            new_parent_time += increment;
            double lb = 0;
            for (int i = 0; i <= std::max(sibling, detach) - n; i++) {
                lb += t[i];
            }
            if (parent == n - 2) {
                log_alpha -= old_parent_time - lb;
            } else {
                double old_ub = t[0] + t[1];
                int ind = 1;
                while (left_child[ind] != n + parent &&
                       right_child[ind] != n + parent) {
                    ind++;
                    old_ub += t[ind];
                }
                log_alpha -= log(old_ub - lb);
            }
            sibling = detach_reattach(detach, attach_above, new_parent_time);
            log_alpha += log_likelihood_times() + log_likelihood_types();
            if (log(gsl_rng_uniform_pos(gen)) > log_alpha) {
                ret = 0;
                sibling = detach_reattach(detach, sibling, old_parent_time);
            }
        }
        return ret;
    }

    int metropolis_mutation_rate() {
        int ret = 1;
        double old_rate = mutation_rate;
        double log_alpha = -log_likelihood_types();
        mutation_rate += gsl_ran_gaussian(gen, hybrid_mutation_sd);
        if (mutation_rate < 0) {
            mutation_rate = -mutation_rate;
        }
        log_alpha += log_likelihood_types();
        if (log(gsl_rng_uniform_pos(gen)) > log_alpha) {
            ret = 0;
            mutation_rate = old_rate;
        }
        return ret;
    }

    void pivot_left(const int i) {
        int tmp = 0;
        tmp = right_child[i - 1];
        if (left_child[i] == n + i - 1) {
            right_child[i - 1] = right_child[i];
            right_child[i] = tmp;
        } else {
            right_child[i - 1] = left_child[i];
            left_child[i] = tmp;
        }
        return;
    }

    void pivot_right(const int i) {
        int tmp = 0;
        tmp = left_child[i - 1];
        if (left_child[i] == n + i - 1) {
            left_child[i - 1] = right_child[i];
            right_child[i] = tmp;
        } else {
            left_child[i - 1] = left_child[i];
            left_child[i] = tmp;
        }
        return;
    }

    void swap(const int i) {
        for (int j = i + 1; j < n - 1; j++) {
            if (left_child[j] == n + i - 1) {
                left_child[j]++;
            } else if (left_child[j] == n + i) {
                left_child[j]--;
            }
            if (right_child[j] == n + i - 1) {
                right_child[j]++;
            } else if (right_child[j] == n + i) {
                right_child[j]--;
            }
        }
        int tmp = left_child[i];
        left_child[i] = left_child[i - 1];
        left_child[i - 1] = tmp;
        tmp = right_child[i];
        right_child[i] = right_child[i - 1];
        right_child[i - 1] = tmp;
        return;
    }

    void run(const double t_max, double &hybrid_spr_acceptance_rate,
             double &hybrid_mutation_acceptance_rate) {
        double sim_time = 0;
        double time_increment, tmp_increment, c;
        double max_drop = 5;
        int ind = 0;
        int metropolis_steps = 0;
        std::vector<int> branches;
        while (sim_time < t_max) {
            branches.resize(n);
            for (int i = 0; i < n; i++) {
                branches[i] = i;
            }
            if (hybrid_rate > 0) {
                time_increment = gsl_ran_exponential(gen, 1 / hybrid_rate);
                ind = -2; // indicates a metropolis step
            } else {
                // there must be a finite maximum step size in the rare case
                // that all velocities are positive.
                time_increment = 2 * mutation_rate;
                // ind = -1 indicates that the time truncation has been hit
                // without a velocity flip or boundary crossing.
                ind = -1;
            }
            if (v[0] < 0 && -t[0] / v[0] < time_increment) {
                if (leaf_to_row(left_child[0]) == leaf_to_row(right_child[0])) {
                    c = 1;
                    ind = 0;
                } else {
                    c = max_drop;
                }
                time_increment = -t[0] / (c * v[0]);
            }
            for (int i = 1; i < n - 1; i++) {
                if (v[i] < 0 && -t[i] / v[i] < time_increment) {
                    time_increment = -t[i] / v[i];
                    ind = i;
                }
            }
            if (mutation_vel < 0 &&
                -mutation_rate / (max_drop * mutation_vel) < time_increment) {
                time_increment = -mutation_rate / (max_drop * mutation_vel);
                ind = -1;
            }
            fill_bounds(time_increment);
            for (int i = 0; i < n; i++) {
                tmp_increment = sample_increment(i, time_increment, branches);
                if (tmp_increment < time_increment) {
                    time_increment = tmp_increment;
                    ind = i;
                }
                if (i < n - 1) {
                    branches.erase(std::lower_bound(
                        branches.begin(), branches.end(), left_child[i]));
                    branches.erase(std::lower_bound(
                        branches.begin(), branches.end(), right_child[i]));
                    branches.push_back(n + i);
                }
            }
            for (int i = 0; i < n - 1; i++) {
                t[i] += v[i] * time_increment;
            }
            mutation_rate += mutation_vel * time_increment;
            sim_time += time_increment;
            if (ind > 0 && ind < n - 1 && t[ind] <= 0) {
                t[ind] = 0;
                if (left_child[ind] == n + ind - 1 ||
                    right_child[ind] == n + ind - 1) {
                    if (gsl_rng_uniform(gen) < 0.5) {
                        pivot_right(ind);
                    } else {
                        pivot_left(ind);
                    }
                } else {
                    swap(ind);
                }
            }
            if (ind >= 0 && ind < n - 1) {
                v[ind] *= -1;
            } else if (ind == n - 1) {
                mutation_vel *= -1;
            } else if (ind == -2) {
                std::cout << mutation_rate * nsites << " " << tree_height()
                          << " 0" << std::endl;
                metropolis_steps++;
                hybrid_spr_acceptance_rate += metropolis_spr();
                hybrid_mutation_acceptance_rate += metropolis_mutation_rate();
            }
            std::cout << mutation_rate * nsites << " " << tree_height();
            if (hybrid_rate > 0) {
                if (ind == -2) {
                    std::cout << " 1";
                } else {
                    std::cout << " 0";
                }
            }
            std::cout << std::endl;
        }
        hybrid_spr_acceptance_rate /= metropolis_steps;
        hybrid_mutation_acceptance_rate /= metropolis_steps;
        return;
    }

    double tree_height() const {
        double ret = 0;
        for (int i = 0; i < n - 1; i++) {
            ret += t[i];
        }
        return ret;
    }

    int n, nsites;
    double hybrid_rate, mutation_rate, mutation_vel, hybrid_mutation_sd;
    std::vector<int> row_count, left_child, right_child;
    std::vector<double> t, v;
    std::vector<std::vector<std::vector<double>>> l, L, m, M, L_tmp, M_tmp;
    gsl_spmatrix *data;
    gsl_rng *gen;
};

#endif
