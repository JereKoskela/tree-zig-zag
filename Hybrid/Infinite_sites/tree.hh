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
#include <list>
#include <sstream>
#include <unistd.h>
#include <vector>

struct Sim {

    Sim(const char *filename)
    : n(0), v_theta(), theta(), t_max(), sd_mut(), mh_rate(), row_count(), col_sum(), left_child(), right_child(), nmut(), t(), v()
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
        cfg.lookupValue("sd_mutation_rate", sd_mut);
        cfg.lookupValue("mh_rate", mh_rate);
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
                                active.insert(active.begin() + ub_j, active[lb_i]);
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

    double log_likelihood_times() const {
        double ret = 0;
        for (int i = 0; i < n - 1; i++) {
            ret -= gsl_sf_choose(n - i, 2) * t[i];
        }
        return ret;
    }

    double log_likelihood_types() const {
        double ret = 0;
        double bl = 0;
        for (int i = 0; i < 2 * n - 2; i++) {
            bl = branch_length(i, 0);
            if (nmut[i] > 0) {
                ret += nmut[i] * log(theta * bl / 2) - theta * bl / 2 - gsl_sf_lnfact(nmut[i]);
            } else {
                ret -= theta * bl / 2;
            }
        }
        return ret;
    }

    int mutations_in_between(int i, int j) const {
        int ret = 0;
        if (j < n) {
            j = 2 * n - 2;
        }
        int ind = std::min(i, j) + 1 - n;
        while (i != j) {
            if (left_child[ind] == i || right_child[ind] == i) {
                ret += nmut[i];
                i = n + ind;
            }
            if (left_child[ind] == j || right_child[ind] == j) {
                ret += nmut[j];
                j = n + ind;
            }
            assert(ind < n - 1);
            ind++;
        }
        return ret;
    }

    int detach_reattach(const int detach, int attach_above, const double new_time) {
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
        nmut.erase(nmut.begin() + n + parent);
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
        nmut.insert(nmut.begin() + n + ind, number_of_mutations(n + ind));
        nmut[sibling] = number_of_mutations(sibling);
        nmut[attach_above] = number_of_mutations(attach_above);
        return sibling;
    }

    int spr() {
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
        } else if (mutations_in_between(n + parent, n + new_parent) > 0) {
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
                for (int i = std::max(std::max(attach_above, detach) - n, -1) + 1; i <= new_parent; i++) {
                    ub += t[i];
                }
                increment = gsl_ran_flat(gen, new_parent_time, ub) - new_parent_time;
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
                while (left_child[ind] != n + parent && right_child[ind] != n + parent) {
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

    int mh_theta() {
        int ret = 1;
        double old_theta = theta;
        double log_alpha = -log_likelihood_types();
        theta += gsl_ran_gaussian(gen, sd_mut);
        if (theta < 0) {
            theta = -theta;
        }
        log_alpha += log_likelihood_types();
        if (log(gsl_rng_uniform_pos(gen)) > log_alpha) {
            theta = old_theta;
            ret = 0;
        }
        return ret;
    }

    double sample_increment(const int i, const double t_max, std::vector<int> &branches) {
        double ret = 0;
        double rate = 0;
        if (i == -1) {
            rate = -(double)col_sum.size() / (theta + v_theta * t_max);
            if (v_theta > 0) {
                for (int j = 0; j < n - 1; j++) {
                    rate += (n - j) * (t[j] + fmax(v[j] * t_max, 0)) / 2;
                }
            } else {
                for (int j = 0; j < n - 1; j++) {
                    rate += (n - j) * (t[j] + fmin(v[j] * t_max, 0)) / 2;
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
                rate = (n - i) * (n + theta + fmax(v_theta * t_max, 0) - i - 1) / 2;
                for (unsigned int j = 0; j < branches.size(); j++) {
                    if (nmut[branches[j]] > 0 && branch_velocity(branches[j]) > 0) {
                        rate -= nmut[branches[j]] / branch_length(branches[j], t_max);
                    } else if (nmut[branches[j]] > 0) {
                        rate -= nmut[branches[j]] / branch_length(branches[j], 0);
                    }
                }
            } else {
                rate = (n - i) * (n + theta + fmin(v_theta * t_max, 0) - i - 1) / 2;
                for (unsigned int j = 0; j < branches.size(); j++) {
                    if (nmut[branches[j]] > 0 && branch_velocity(branches[j]) > 0) {
                        rate -= nmut[branches[j]] / branch_length(branches[j], 0);
                    } else if (nmut[branches[j]] > 0) {
                        rate -= nmut[branches[j]] / branch_length(branches[j], t_max);
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
                ret = t_max;
            }
            if (ret >= t_max) {
                ret = t_max;
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
        int tmp = right_child[i - 1];
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
        int tmp = left_child[i - 1];
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

    void run(double &accept_spr, double &accept_theta) {
        double sim_time = 0;
        double time_increment, tmp_increment, c;
        int ind = 0;
        int mh_steps = 0;
        int ac_spr = 0;
        int ac_theta = 0;
        std::vector<int> branches;
        while (sim_time < t_max) {
            time_increment = gsl_ran_exponential(gen, 1 / mh_rate);
            ind = -3;
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
                ind = -1;
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
            } else if (ind == -3) {
                std::cout << theta << " " << tree_height() << " " << total_branch_length() << " 1" << std::endl;
                mh_steps++;
                ac_spr += spr();
                ac_theta += mh_theta();
            }
            std::cout << theta << " " << tree_height() << " " << total_branch_length() << " 0" << std::endl;
        }
        accept_spr = (double)ac_spr / mh_steps;
        accept_theta = (double)ac_theta / mh_steps;
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
    double v_theta, theta, t_max, sd_mut, mh_rate;
    std::vector<int> row_count, col_sum, left_child, right_child, nmut;
    std::vector<double> t, v;
    gsl_spmatrix *data;
    gsl_rng *gen;
};

#endif
