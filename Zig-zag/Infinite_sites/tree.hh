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
    : n(0), v_theta(), theta(), t_max(), localised(), row_count(), col_sum(), left_child(), right_child(), nmut(), t(), v()
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
        cfg.lookupValue("localised", localised);
        std::string line, token;
        int l = 0;
        std::ifstream file;
        file.open(datafile);
        while(getline(file, line)) {
            std::vector<int> row;
            std::stringstream iss;
            iss << line;
            while (getline(iss, token, ' ')) {
                row.push_back(atoi(token.c_str()));
            }
            if (l == 0) {
                col_sum.resize(row.size() - 1, 0);
            }
            for (unsigned int i = 0; i < row.size() - 1; i++) {
                if (row[i] == 1) {
                    gsl_spmatrix_set(data, l, i, 1);
                    col_sum[i] += row.back();
                }
            }
            row_count.push_back(row.back());
            n += row.back();
            l++;
        }
        left_child.resize(n - 1, -1);
        right_child.resize(n - 1, -1);
        v.resize(n - 1, 1);
        for (int i = 0; i < n - 1; i++) {
            v[i] = 1 / gsl_sf_choose(n - i, 2);
            if (gsl_rng_uniform(gen) < 0.5) {
                v[i] *= -1;
            }
        }
        if (gsl_rng_uniform(gen) < 0.5) {
            v_theta *= -1;
        }
        t.resize(n - 1, 0);
        nmut.resize(2 * n - 1, -1);
        theta = watterson_estimator();
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
        return col_sum.size() / harmonic;
    }

    void simplify(gsl_matrix *tmp_data, std::vector<int> &tmp_row_counts,
                  std::vector<int> &tmp_col_sums, std::vector<int> &active) {
        bool changed = false;
        int diff = 0, ub_j = 0, lb_i = 0;
        for (unsigned int i = 0; i < tmp_row_counts.size(); i++) {
            changed = false;
            if (tmp_row_counts[i] == 1) {
                for (unsigned int j = 0; j < tmp_col_sums.size(); j++) {
                    if (gsl_matrix_get(tmp_data, i, j) == 1 && tmp_col_sums[j] == 1) {
                        gsl_matrix_set(tmp_data, i, j, 0);
                        tmp_col_sums[j]--;
                        changed = true;
                    }
                }
            }
            if (changed) {
                for (unsigned int j = 0; j < tmp_row_counts.size(); j++) {
                    if (j != i && tmp_row_counts[j] > 0) {
                        diff = 0;
                        for (unsigned int k = 0; k < tmp_col_sums.size(); k++) {
                            diff += std::abs(gsl_matrix_get(tmp_data, i, k) - gsl_matrix_get(tmp_data, j, k));
                        }
                        if (diff == 0) {
                            if (std::abs(int(j) - int(i)) > 1) {
                                ub_j = 0;
                                for (unsigned int k = 0; k <= j; k++) {
                                    ub_j += tmp_row_counts[k];
                                }
                                lb_i = 0;
                                for (unsigned int k = 0; k < i; k++) {
                                    lb_i += tmp_row_counts[k];
                                }
                                active.insert(active.begin() + ub_j,
                                              active[lb_i]);
                                if (j < i) {
                                    active.erase(active.begin() + lb_i + 1, active.begin() + lb_i + 2);
                                } else {
                                    active.erase(active.begin() + lb_i, active.begin() + lb_i + 1);
                                }
                            }
                            tmp_row_counts[j] += tmp_row_counts[i];
                            tmp_row_counts[i] = 0;
                            for (unsigned int k = 0; k < tmp_col_sums.size(); k++) {
                                gsl_matrix_set(tmp_data, i, k, 0);
                            }
                            break;
                        }
                    }
                }
            }
        }
        return;
    }

    void sample_tree() {
        gsl_matrix *tmp_data = gsl_matrix_alloc(data->size1, data->size2);
        gsl_spmatrix_sp2d(tmp_data, data);
        std::vector<int> tmp_row_count = row_count, tmp_col_sum = col_sum;
        std::vector<int> active(n);
        for (int i = 0; i < n; i++) {
            active[i] = i;
        }
        simplify(tmp_data, tmp_row_count, tmp_col_sum, active);
        int next_parent = n;
        std::vector<double> p(tmp_row_count.size(), 0);
        double sum_p = 0, ub = 0, coin = 0;
        int least_index = 0, c1 = 0, c2 = 0, row = 0;
        while (active.size() > 1) {
            sum_p = 0;
            for (unsigned int i = 0; i < tmp_row_count.size(); i++) {
                if (tmp_row_count[i] > 1) {
                    p[i] = gsl_sf_choose(tmp_row_count[i], 2);
                } else {
                    p[i] = 0;
                }
                sum_p += p[i];
            }
            row = 0;
            ub = p[0];
            coin = gsl_rng_uniform(gen) * sum_p;
            while (coin > ub) {
                row++;
                ub += p[row];
            }
            least_index = 0;
            for (int i = 0; i < row; i++) {
                least_index += tmp_row_count[i];
            }
            c1 = least_index + gsl_rng_uniform_int(gen, tmp_row_count[row]);
            do {
                c2 = least_index + gsl_rng_uniform_int(gen, tmp_row_count[row]);
            } while (c1 == c2);
            left_child[next_parent - n] = active[c1];
            right_child[next_parent - n] = active[c2];
            t[next_parent - n] = gsl_ran_exponential(gen, 1 / gsl_sf_choose(active.size(), 2));
            tmp_row_count[row]--;
            active[std::min(c1, c2)] = next_parent;
            active.erase(active.begin() + std::max(c1, c2));
            next_parent++;
            for (unsigned int j = 0; j < tmp_col_sum.size(); j++) {
                tmp_col_sum[j] -= gsl_matrix_get(tmp_data, row, j);
            }
            if (tmp_row_count[row] == 1) {
                simplify(tmp_data, tmp_row_count, tmp_col_sum, active);
            }
        }
        gsl_matrix_free(tmp_data);
        for (int i = 0; i < 2 * n - 1; i++) {
            nmut[i] = number_of_mutations(i);
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

    void leaves(std::vector<int> &l, const int i) const {
        if (i < n) {
            l.push_back(i);
        } else {
            leaves(l, left_child[i - n]);
            leaves(l, right_child[i - n]);
        }
        return;
    }

    int number_of_mutations(const int i) const {
        int ret = 0, leaf_sum;
        std::vector<int> l;
        leaves(l, i);
        for (unsigned int site = 0; site < col_sum.size(); site++) {
            if ((int)l.size() == col_sum[site]) {
                leaf_sum = 0;
                for (unsigned int j = 0; j < l.size(); j++) {
                    if (gsl_spmatrix_get(data, leaf_to_row(l[j]), site) == 1) {
                        leaf_sum++;
                    }
                }
                if (leaf_sum == col_sum[site]) {
                    ret++;
                }
            }
        }
        return ret;
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

    double sample_t_up(const int i, std::vector<int> &branches) const {
        double ret = 0;
        double tmp, discriminant, root_1, root_2, alpha = 0, u;
        if ((int)branches.size() > n - i) {
            branches.erase(std::lower_bound(branches.begin(), branches.end(), left_child[i - 1]));
            branches.erase(std::lower_bound(branches.begin(), branches.end(), right_child[i - 1]));
            branches.push_back(n + i - 1);
        }
        do {
            u = gsl_rng_uniform(gen);
            tmp = (n + theta - i - 1 + v_theta * ret);
            discriminant = tmp * tmp - 4 * v_theta * log(u) / (v[i] * (n - i));
            if (discriminant < 0) {
                ret = std::numeric_limits<double>::max();
                alpha = 1;
            } else {
                root_1 = (sqrt(discriminant) - tmp) / v_theta;
                root_2 = -(sqrt(discriminant) + tmp) / v_theta;
                if (root_1 > 0 && root_2 > 0) {
                    ret += fmin(root_1, root_2);
                } else {
                    ret += fmax(root_1, root_2);
                }
                for (int j = 0; j < n - 1; j++) {
                    if (t[j] + v[j] * ret <= 1e-10) {
                        alpha = 1;
                    }
                }
                if (theta + v_theta * ret <= 1e-10) {
                    alpha = 1;
                }
                if (alpha < 1) {
                    alpha = (n - i) * (n + theta + v_theta * ret - i - 1) / 2;
                    tmp = 0;
                    for (unsigned int j = 0; j < branches.size(); j++) {
                        tmp += nmut[branches[j]] / branch_length(branches[j], ret);
                    }
                    alpha = (alpha - tmp) / alpha;
                }
            }
        } while (gsl_rng_uniform(gen) > alpha);
        return ret;
    }

    double sample_t_down(const int i, std::vector<int> &branches) const {
        double ret = 0;
        double alpha = 0, tmp, u, vel, min_increment, bl;
        int m;
        if ((int)branches.size() > n - i) {
            branches.erase(std::lower_bound(branches.begin(), branches.end(), left_child[i - 1]));
            branches.erase(std::lower_bound(branches.begin(), branches.end(), right_child[i - 1]));
            branches.push_back(n + i - 1);
        }
        do {
            min_increment = std::numeric_limits<double>::max();
            for (unsigned int j = 0; j < branches.size(); j++) {
                u = gsl_rng_uniform_pos(gen);
                m = nmut[branches[j]];
                vel = branch_velocity(branches[j]);
                bl = branch_length(branches[j], ret);
                if (m == 0) {
                    tmp = -t[i] / v[i] - ret;
                } else if (fabs(vel) < 1e-10) {
                    tmp = bl * log(u) / (m * v[i]);
                } else {
                    tmp = bl * (pow(u, vel / (v[i] * m)) - 1) / vel;
                }
                if (tmp < min_increment) {
                    min_increment = tmp;
                }
            }
            ret += min_increment;
            for (int j = 0; j < n - 1; j++) {
                if (t[j] + v[j] * ret <= 1e-10) {
                    alpha = 1;
                    if (j == i) {
                        ret = -t[i] / v[i];
                    }
                }
            }
            if (theta + v_theta * ret <= 1e-10) {
                alpha = 1;
            }
            if (alpha < 1) {
                alpha = 0;
                for (unsigned int j = 0; j < branches.size(); j++) {
                    alpha += nmut[branches[j]] / branch_length(branches[j], ret);
                }
                alpha = (alpha - (n - i) * (n + theta + v_theta * ret - i - 1) / 2) / alpha;
            }
        } while (gsl_rng_uniform(gen) > alpha);
        return ret;
    }

    double sample_theta_up() const {
        double ret = 0;
        double alpha = 0, u, tmp, root_1, root_2, discriminant;
        double sum_v = 0;
        for (int i = 0; i < n - 1; i++) {
            sum_v += (n - i) * v[i];
        }
        do {
            u = gsl_rng_uniform(gen);
            tmp = 0;
            for (int i = 0; i < n - 1; i++) {
                tmp += (n - i) * (t[i] + v[i] * ret);
            }
            discriminant = tmp * tmp - 4 * sum_v * log(u) / v_theta;
            if (discriminant < 0) {
                ret = std::numeric_limits<double>::max();
                alpha = 1;
            } else {
                root_1 = (-tmp + sqrt(discriminant)) / sum_v;
                root_2 = -(tmp + sqrt(discriminant)) / sum_v;
                if (root_1 > 0 && root_2 > 0) {
                    ret += fmin(root_1, root_2);
                } else {
                    ret += fmax(root_1, root_2);
                }
            }
            for (int i = 0; i < n - 1; i++) {
                if (t[i] + v[i] * ret <= 1e-10) {
                    alpha = 1;
                }
            }
            if (alpha < 1) {
                alpha = 0;
                for (int i = 0; i < n - 1; i++) {
                    alpha += (n - i) * (t[i] + v[i] * ret) / 2;
                }
                alpha = (alpha - col_sum.size() / (theta + v_theta * ret)) / alpha;
            }
        } while (gsl_rng_uniform(gen) > alpha);
        return ret;
    }

    double sample_theta_down() const {
        double ret = 0;
        double alpha = 0, u;
        do {
            u = gsl_rng_uniform(gen);
            if (col_sum.size() == 0) {
                ret = -theta / v_theta;
                alpha = 1;
            } else {
                ret += (theta + v_theta * ret) * (pow(u, 1.0 / col_sum.size()) - 1) / v_theta;
            }
            for (int i = 0; i < n - 1; i++) {
                if (t[i] + v[i] * ret <= 1e-10) {
                    alpha = 1;
                }
            }
            if (theta + v_theta * ret <= 1e-10) {
                ret = -theta / v_theta;
                alpha = 1;
            }
            if (alpha < 1) {
                alpha = 0;
                for (int i = 0; i < n - 1; i++) {
                    alpha += (n - i) * (t[i] + v[i] * ret) / 2;
                }
                alpha = (col_sum.size() / (theta + v_theta * ret) - alpha) / (col_sum.size() / (theta + v_theta * ret));
            }
        } while (gsl_rng_uniform(gen) > alpha);
        return ret;
    }

    double sample_localised_increment(const int i, const double ub, std::vector<int> &branches) {
        double ret = 0;
        double rate = 0;
        if (i == -1) {
            rate = -(double)col_sum.size() / (theta + v_theta * ub);
            if (v_theta > 0) {
                for (int j = 0; j < n - 1; j++) {
                    rate += (n - j) * (t[j] + fmax(v[j] * ub, 0)) / 2;
                }
            } else {
                for (int j = 0; j < n - 1; j++) {
                    rate += (n - j) * (t[j] + fmin(v[j] * ub, 0)) / 2;
                }
            }
            rate *= v_theta;
        } else {
            if ((int)branches.size() > n - i) {
                branches.erase(std::lower_bound(branches.begin(), branches.end(), left_child[i - 1]));
                branches.erase(std::lower_bound(branches.begin(), branches.end(), right_child[i - 1]));
                branches.push_back(n + i - 1);
            }
            if (v[i] > 0) {
                rate = (n - i) * (n + theta + fmax(v_theta * ub, 0) - i - 1) / 2;
                for (unsigned int j = 0; j < branches.size(); j++) {
                    if (nmut[branches[j]] > 0 && branch_velocity(branches[j]) > 0) {
                        rate -= nmut[branches[j]] / branch_length(branches[j], ub);
                    } else if (nmut[branches[j]] > 0) {
                        rate -= nmut[branches[j]] / branch_length(branches[j], 0);
                    }
                }
            } else {
                rate = (n - i) * (n + theta + fmin(v_theta * ub, 0) - i - 1) / 2;
                for (unsigned int j = 0; j < branches.size(); j++) {
                    if (nmut[branches[j]] > 0 && branch_velocity(branches[j]) > 0) {
                        rate -= nmut[branches[j]] / branch_length(branches[j], 0);
                    } else if (nmut[branches[j]] > 0) {
                        rate -= nmut[branches[j]] / branch_length(branches[j], ub);
                    }
                }
            }
            rate *= v[i];
        }
        double alpha;
        do {
            if (rate > 0) {
                ret += gsl_ran_exponential(gen, 1 / rate);
            } else {
                ret = ub;
            }
            if (ret >= ub) {
                ret = ub;
                alpha = 1;
            } else {
                if (i == -1) {
                    alpha = -(double)col_sum.size() / (theta + v_theta * ret);
                    for (int j = 0; j < n - 1; j++) {
                        alpha += (n - j) * (t[j] + v[j] * ret) / 2;
                    }
                    alpha *= v_theta / rate;
                } else {
                    alpha = (n - i) * (n + theta + v_theta * ret - i - 1) / 2;
                    for (unsigned int j = 0; j < branches.size(); j++) {
                        if (nmut[branches[j]] > 0) {
                            alpha -= nmut[branches[j]] / branch_length(branches[j], ret);
                        }
                    }
                    alpha *= v[i] / rate;
                }
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
        tmp = nmut[n + i];
        nmut[n + i] = nmut[n + i - 1];
        nmut[n + i - 1] = tmp;
        return;
    }

    void run() {
        double sim_time = 0;
        double time_increment, tmp_increment;
        int ind = 0;
        int c = 0;
        std::vector<int> branches;
        while (sim_time < t_max) {
            time_increment = std::numeric_limits<double>::max();
            branches.resize(n);
            for (int i = 0; i < n; i++) {
                branches[i] = i;
            }
            if (localised) {
                time_increment = 2 * theta;
                ind = -2; // records arrival at the end of the time window - no velocity flip needed
                for (int i = 0; i < n - 1; i++) {
                    if (v[i] < 0) {
                        c = 1;
                        if (i == 0) {
                            if (nmut[left_child[0]] > 0 || nmut[right_child[0]] > 0) {
                                c = 5;
                            }
                        } else {
                            if ((n + i - 1 == left_child[i] || n + i - 1 == right_child[i]) && nmut[n + i - 1] > 0) {
                                c = 5;
                            }
                        }
                        if (-t[i] / (c * v[i]) < time_increment) {
                            time_increment = -t[i] / (c * v[i]);
                            if (c == 1) {
                                ind = i;
                            } else {
                                ind = -2;
                            }
                        }
                    }
                }
                if (v_theta < 0 && -theta / (5 * v_theta) < time_increment) {
                    time_increment = -theta / (5 * v_theta);
                    ind = -2;
                }
                for (int i = 0; i < n - 1; i++) {
                    tmp_increment = sample_localised_increment(i, time_increment, branches);
                    if (tmp_increment < time_increment) {
                        time_increment = tmp_increment;
                        ind = i;
                    }
                }
                tmp_increment = sample_localised_increment(-1, time_increment, branches);
                if (tmp_increment < time_increment) {
                    time_increment = tmp_increment;
                    ind = -1;
                }
            } else {
                for (int i = 0; i < n - 1; i++) {
                    if (v[i] > 0) {
                        tmp_increment = sample_t_up(i, branches);
                    } else {
                        tmp_increment = sample_t_down(i, branches);
                    }
                    if (tmp_increment < time_increment) {
                        time_increment = tmp_increment;
                        ind = i;
                    }
                }
                if (v_theta > 0) {
                    tmp_increment = sample_theta_up();
                } else {
                    tmp_increment = sample_theta_down();
                }
                if (tmp_increment < time_increment) {
                    time_increment = tmp_increment;
                    ind = -1; // records that the next flip is in the theta direction
                }
            }
            for (int i = 0; i < n - 1; i++) {
                t[i] += v[i] * time_increment;
            }
            theta += v_theta * time_increment;
            sim_time += time_increment;
            if (ind > 0 && t[ind] <= 0) {
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
            assert((ind == -2 && localised) || ind > -2);
            std::cout << theta << " " << tree_height() << " " << total_branch_length() << std::endl;
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

    int n;
    double v_theta, theta, t_max;
    bool localised;
    std::vector<int> row_count, col_sum, left_child, right_child, nmut;
    std::vector<double> t, v;
    gsl_spmatrix *data;
    gsl_rng *gen;
};

#endif
