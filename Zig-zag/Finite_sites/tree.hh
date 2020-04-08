#ifndef SIM
#define SIM

#include <algorithm>
#include <cassert>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_spmatrix.h>
#include <libconfig.h++>
#include <limits>
#include <sstream>
#include <unistd.h>
#include <vector>

struct Sim {

    Sim(const char *filename)
    : n(0), theta(), v_theta(), t_max(), row_count(), left_child(), right_child(), t(), v(), l(), L(), m(), M(), L_tmp(), M_tmp()
    {
        gen = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(gen, time(NULL) * getpid());
        data = gsl_spmatrix_alloc(1, 1);
        libconfig::Config cfg;
        cfg.readFile(filename);
        std::string datafile;
        cfg.lookupValue("datafile", datafile);
        cfg.lookupValue("v_theta", v_theta);
        cfg.lookupValue("run_length", t_max);
        std::string line, token;
        int s = 0;
        std::ifstream file;
        file.open(datafile);
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
        v_theta /= nsites;
        if (gsl_rng_uniform(gen) < 0.5) {
            v_theta *= -1;
        }
        t.resize(n - 1, 0);
        std::vector<double> tmp(2);
        std::vector<std::vector<double> > tmp_2(2 * n - 1, tmp);
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
        theta = watterson_estimator() / nsites;
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
                t[next_parent - n] = gsl_ran_exponential(gen, 1 / gsl_sf_choose(active.size(), 2));
                active[i] = next_parent;
                active.erase(active.begin() + i + 1);
                next_parent++;
            }
        }
        for (int i = 0; i < (int)row_count.size() - 1; i++) {
            left_child[next_parent - n] = active[0];
            right_child[next_parent - n] = active[1];
            t[next_parent - n] = gsl_ran_exponential(gen, 1 / gsl_sf_choose(active.size(), 2));
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

    double post_order(const int i, const int c, const int h, const double s) const {
        double bl = branch_length(c, s);
        int g = (h + 1) % 2;
        double ret = (L_tmp[i][c][0] + L_tmp[i][c][1]
            + (L_tmp[i][c][h] - L_tmp[i][c][g]) * exp(-(theta + v_theta * s) * bl)) / 2;
        return ret;
    }

    double post_order_ub(const int i, const int c, const int h, const double s) const {
        double ret = L[i][c][0] + L[i][c][1];
        double bl = branch_length(c, 0);
        double vl = branch_velocity(c);
        int g = (h + 1) % 2;
        if (L[i][c][h] > L[i][c][g]) {
            ret += (L[i][c][h] - L[i][c][g]) * exp(-fmin(theta * bl, (theta + v_theta * s) * (bl + vl * s)));
        } else {
            double crit = -bl / (2 * vl) - theta / (2 * v_theta);
            if (crit > 0 && crit < s) {
                ret += (L[i][c][h] - L[i][c][g]) * exp(-(theta - v_theta * bl / vl) * (bl - vl * theta / v_theta) / 4);
            } else {
                ret += (L[i][c][h] - L[i][c][g]) * exp(-fmax(theta * bl, (theta + v_theta * s) * (bl + vl * s)));
            }
        }
        return fmin(ret / 2, 1);
    }

    double post_order_lb(const int i, const int c, const int h, const double s) const {
        double ret = l[i][c][0] + l[i][c][1];
        double bl = branch_length(c, 0);
        double vl = branch_velocity(c);
        int g = (h + 1) % 2;
        if (l[i][c][h] > l[i][c][g]) {
            double crit = -bl / (2 * vl) - theta / (2 * v_theta);
            if (crit > 0 && crit < s) {
                ret += (l[i][c][h] - l[i][c][g]) * exp(-(theta - v_theta * bl / vl) * (bl - vl * theta / v_theta) / 4);
            } else {
                ret += (l[i][c][h] - l[i][c][g]) * exp(-fmax(theta * bl, (theta + v_theta * s) * (bl + vl * s)));
            }
        } else {
            ret += (l[i][c][h] - l[i][c][g]) * exp(-fmin(theta * bl, (theta + v_theta * s) * (bl + vl * s)));
        }
        return fmax(ret / 2, 0);
    }

    double pre_order(const int i, const int p, const int c, const int d, const int h, const double s) const {
        double bl = branch_length(c, s);
        int g = (h + 1) % 2;
        double c0 = M_tmp[i][p][h] * post_order(i, d, h, s);
        double c1 = M_tmp[i][p][g] * post_order(i, d, g, s);
        double ret = (c0 + c1 + (c0 - c1) * exp(-(theta + v_theta * s) * bl)) / 2;
        return ret;
    }

    double pre_order_ub(const int i, const int p, const int c, const int d, const int h, const double s) const {
        double bl = branch_length(c, 0);
        double vl = branch_velocity(c);
        int g = (h + 1) % 2;
        double c0 = M[i][p][h] * post_order_ub(i, d, h, s);
        double c1 = M[i][p][g] * post_order_ub(i, d, g, s);
        double ret = c0 + c1;
        if (c0 > c1) {
            ret += (c0 - c1) * exp(-fmin(theta * bl, (theta + v_theta * s) * (bl + vl * s)));
        } else {
            double crit = -bl / (2 * vl) - theta / (2 * v_theta);
            if (crit > 0 && crit < s) {
                ret += (c0 - c1) * exp(-(theta - v_theta * bl / vl) * (bl - vl * theta / v_theta) / 4);
            } else {
                ret += (c0 - c1) * exp(-fmax(theta * bl, (theta + v_theta * s) * (bl + vl * s)));
            }
        }
        return fmin(ret / 2, 1);
    }

    double pre_order_lb(const int i, const int p, const int c, const int d, const int h, const double s) const {
        double bl = branch_length(c, 0);
        double vl = branch_velocity(c);
        int g = (h + 1) % 2;
        double c0 = m[i][p][h] * post_order_lb(i, d, h, s);
        double c1 = m[i][p][g] * post_order_lb(i, d, g, s);
        double ret = c0 + c1;
        if (c0 > c1) {
            double crit = -bl / (2 * vl) - theta / (2 * v_theta);
            if (crit > 0 && crit < s) {
                ret += (c0 - c1) * exp(-(theta - v_theta * bl / vl) * (bl - vl * theta / v_theta) / 4);
            } else {
                ret += (c0 - c1) * exp(-fmax(theta * bl, (theta + v_theta * s) * (bl + vl * s)));
            }
        } else {
            ret += (c0 - c1) * exp(-fmin(theta * bl, (theta + v_theta * s) * (bl + vl * s)));
        }
        return fmax(ret / 2, 0);
    }

    void fill_bounds(const double s) {
        int c, d;
        for (int i = 0; i < nsites; i++) {
            for (int j = n; j < 2 * n - 1; j++) {
                L[i][j][0] = post_order_ub(i, left_child[j - n], 0, s)
                    * post_order_ub(i, right_child[j - n], 0, s);
                L[i][j][1] = post_order_ub(i, left_child[j - n], 1, s)
                    * post_order_ub(i, right_child[j - n], 1, s);
                l[i][j][0] = post_order_lb(i, left_child[j - n], 0, s)
                    * post_order_lb(i, right_child[j - n], 0, s);
                l[i][j][1] = post_order_lb(i, left_child[j - n], 1, s)
                    * post_order_lb(i, right_child[j - n], 1, s);
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

    double grad_log_L_ub(const int ind, const double s, const std::vector<int> &branches) {
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
                    ret += (theta + fmax(0, v_theta * s)) * temp / (l[i][2 * n - 2][0] + l[i][2 * n - 2][1]);
                } else {
                    ret += (theta + fmin(0, v_theta * s)) * temp / (L[i][2 * n - 2][0] + L[i][2 * n - 2][1]);
                }
            }
        } else {
            double bl, vl, ret_i;
            for (int i = 0; i < nsites; i++) {
                ret_i = 0;
                for (int j = 0; j < 2 * n - 2; j++) {
                    temp = 0;
                    if (v_theta > 0) {
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
                    if (temp * v_theta < 0) {
                        ret_i += (bl + fmax(0, vl * s)) * temp;
                    } else {
                        ret_i += (bl + fmin(0, vl * s)) * temp;
                    }
                }
                if (ret_i * v_theta < 0) {
                    ret += ret_i / (l[i][2 * n - 2][0] + l[i][2 * n - 2][1]);
                } else {
                    ret += ret_i / (L[i][2 * n - 2][0] + L[i][2 * n - 2][1]);
                }
            }
        }
        return ret;
    }

    double grad_log_L(const int ind, const double s, const std::vector<int> &branches) {
        double ret = 0;
        int c, d;
        for (int i = 0; i < nsites; i++) {
            for (int j = n; j < 2 * n - 1; j++) {
                L_tmp[i][j][0] = post_order(i, left_child[j - n], 0, s)
                    * post_order(i, right_child[j - n], 0, s);
                L_tmp[i][j][1] = post_order(i, left_child[j - n], 1, s)
                    * post_order(i, right_child[j - n], 1, s);
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
                    ret += (theta + v_theta * s) * (M_tmp[i][c][0] - M_tmp[i][c][1])
                        * (L_tmp[i][c][1] - L_tmp[i][c][0])
                        / (L_tmp[i][2 * n - 2][0] + L_tmp[i][2 * n - 2][1]);
                }
            } else {
                for (int j = 0; j < 2 * n - 2; j++) {
                    ret += branch_length(j, s) * (M_tmp[i][j][0] - M_tmp[i][j][1])
                        * (L_tmp[i][j][1] - L_tmp[i][j][0])
                        / (L_tmp[i][2 * n - 2][0] + L_tmp[i][2 * n - 2][1]);
                }
            }
        }
        return ret;
    }

    double sample_increment(const int i, const double t_max, std::vector<int> &branches) {
        double ret = 0;
        double rate = 0;
        if (i > -1 && (int)branches.size() > n - i) {
            branches.erase(std::lower_bound(branches.begin(), branches.end(), left_child[i - 1]));
            branches.erase(std::lower_bound(branches.begin(), branches.end(), right_child[i - 1]));
            branches.push_back(n + i - 1);
        }
        if (i == -1) {
            rate = -v_theta * grad_log_L_ub(-1, t_max, branches);
        } else {
            rate = v[i] * (gsl_sf_choose(n - i, 2) - grad_log_L_ub(i, t_max, branches));
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
                if (i == -1) {
                    alpha = -v_theta * grad_log_L(-1, ret, branches) / rate;
                } else {
                    alpha = v[i] * (gsl_sf_choose(n - i, 2) - grad_log_L(i, ret, branches)) / rate;
                }
                assert(alpha < 1);
            }
        } while (gsl_rng_uniform(gen) > alpha);
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

    void run() {
        double sim_time = 0;
        double time_increment, tmp_increment;
        int ind = -2;
        std::vector<int> branches;
        while (sim_time < t_max) {
            //time_increment = std::numeric_limits<double>::max();
            time_increment = 5e-2;
            ind = -2;
            if (v[0] < 0 && -t[0] / v[0] < time_increment) {
                if (leaf_to_row(left_child[0]) == leaf_to_row(right_child[0])) {
                    time_increment = -t[0] / v[0];
                    ind = 0;
                } else {
                    time_increment = -t[0] / (2 * v[0]);
                    ind = -2;
                }
            }
            for (int i = 1; i < n - 1; i++) {
                if (v[i] < 0 && -t[i] / v[i] < time_increment) {
                    time_increment = -t[i] / v[i];
                    ind = i;
                }
            }
            if (v_theta < 0 && -theta / (5 * v_theta) < time_increment) {
                time_increment = -theta / (5 * v_theta);
                ind = -2; // ind = -2 after simulating flip times means the next
                // event is due to theta falling to half of its value, requiring
                // no velocity flips.
            }
            fill_bounds(time_increment);
            branches.resize(n);
            for (int i = 0; i < n; i++) {
                branches[i] = i;
            }
            for (int i = 0; i < n - 1; i++) {
                tmp_increment = sample_increment(i, time_increment, branches);
                if (tmp_increment < time_increment) {
                    time_increment = tmp_increment;
                    ind = i;
                }
            }
            tmp_increment = sample_increment(-1, time_increment, branches);
            if (tmp_increment < time_increment) {
                time_increment = tmp_increment;
                ind = -1; // records that the next flip is in the theta direction
            }
            for (int i = 0; i < n - 1; i++) {
                t[i] += v[i] * time_increment;
            }
            theta += v_theta * time_increment;
            sim_time += time_increment;
            if (ind > 0 && t[ind] <= 0) {
                t[ind] = 0;
                if (left_child[ind] == n + ind - 1 || right_child[ind] == n + ind - 1) {
                    if (gsl_rng_uniform(gen) < 0.5) {
                        pivot_right(ind);
                    } else {
                        pivot_left(ind);
                    }
                } else {
                    swap(ind);
                }
            }
            if (ind >= 0) {
                v[ind] *= -1;
            } else if (ind == -1) {
                v_theta *= -1;
            }
            std::cout << theta * nsites << " " << tree_height() << " " << total_branch_length() << std::endl;
        }
        return;
    }

    double tree_height() const {
        double ret = 0;
        for (int i = 0; i < n - 1; i++) {
            ret += t[i];
        }
        return ret;
    }

    double total_branch_length() const {
        double ret = 0;
        for (int i = 0; i < n - 1; i++) {
            ret += (n - i) * t[i];
        }
        return ret;
    }

    int n, nsites;
    double theta, v_theta, t_max;
    std::vector<int> row_count, left_child, right_child;
    std::vector<double> t, v;
    std::vector<std::vector<std::vector<double> > > l, L, m, M, L_tmp, M_tmp;
    gsl_spmatrix *data;
    gsl_rng *gen;
};

#endif
