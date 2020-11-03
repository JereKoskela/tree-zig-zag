#ifndef SIM
#define SIM

#include <algorithm>
#include <cassert>
#include <fstream>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_spmatrix.h>
#include <iomanip>
#include <sstream>
#include <unistd.h>
#include <vector>

struct Sim {

    Sim(std::string filename)
        : mutation_rate(), ll_times(), ll_types(), sample_size(0), row_counts(),
          col_sums(), parent(), child_1(), child_2(), event_to_parent(), nmut(),
          node_times(), event_times() {
        gen = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(gen, time(NULL) * getpid());
        std::ifstream file;
        file.open(filename);
        std::string line, token;
        gsl_spmatrix *tmp = gsl_spmatrix_alloc(1, 1);
        int l = 0;
        while (getline(file, line)) {
            std::vector<int> tmp_vec;
            std::stringstream iss;
            iss << line;
            while (getline(iss, token, ' ')) {
                tmp_vec.push_back(atoi(token.c_str()));
            }
            if (l == 0) {
                col_sums.resize(tmp_vec.size() - 1, 0);
            }
            for (unsigned int i = 0; i < tmp_vec.size() - 1; i++) {
                if (tmp_vec[i] == 1) {
                    gsl_spmatrix_set(tmp, l, i, 1);
                    col_sums[i] += tmp_vec.back();
                }
            }
            row_counts.push_back(tmp_vec.back());
            sample_size += tmp_vec.back();
            l++;
        }
        parent.resize(2 * sample_size - 1, -1);
        child_1.resize(2 * sample_size - 1, -1);
        child_2.resize(2 * sample_size - 1, -1);
        event_to_parent.resize(sample_size - 1, -1);
        nmut.resize(2 * sample_size - 1, -1);
        node_times.resize(2 * sample_size - 1, 0.0);
        event_times.resize(sample_size - 1, 0.0);
        data = gsl_matrix_alloc(tmp->size1, tmp->size2);
        tmp_data = gsl_matrix_alloc(tmp->size1, tmp->size2);
        gsl_spmatrix_sp2d(data, tmp);
        gsl_spmatrix_free(tmp);
        mutation_rate = watterson_estimator() / 2;
        sample_tree();
        ll_times = log_likelihood_times();
        ll_types = log_likelihood_types();
    }

    Sim(const int n, const double m)
        : mutation_rate(m), ll_times(), ll_types(), sample_size(n),
          row_counts(), col_sums(), parent(2 * n - 1, -1),
          child_1(2 * n - 1, -1), child_2(2 * n - 1, -1),
          event_to_parent(n - 1, -1), node_times(2 * n - 1, 0),
          event_times(n - 1, 0) {
        gen = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(gen, time(NULL) * getpid());
    }

    ~Sim() {
        gsl_rng_free(gen);
        gsl_matrix_free(data);
        gsl_matrix_free(tmp_data);
    }

    double watterson_estimator() const {
        double harmonic = 1;
        for (int i = 2; i < sample_size; i++) {
            harmonic += 1.0 / i;
        }
        return col_sums.size() / harmonic;
    }

    void generate_data() {
        const double theta = mutation_rate;
        int n = 2;
        gsl_spmatrix *tmp = gsl_spmatrix_alloc(1, 1);
        row_counts.push_back(2);
        int child = 0;
        int ub = row_counts[0];
        double coin = 0;
        int next_segregating_site = 0;
        while (n < sample_size + 1) {
            child = 0;
            ub = row_counts[0];
            coin = gsl_rng_uniform(gen);
            while (coin > double(ub) / double(n)) {
                child++;
                ub += row_counts[child];
            }
            if (gsl_rng_uniform(gen) < theta / (n - 1 + theta)) {
                if (row_counts[child] == 1) {
                    gsl_spmatrix_set(tmp, child, next_segregating_site, 1);
                } else {
                    row_counts[child]--;
                    row_counts.push_back(1);
                    int new_child = row_counts.size() - 1;
                    for (int i = 0; i < next_segregating_site; i++) {
                        gsl_spmatrix_set(tmp, new_child, i,
                                         gsl_spmatrix_get(tmp, child, i));
                    }
                    gsl_spmatrix_set(tmp, new_child, next_segregating_site, 1);
                }
                col_sums.push_back(1);
                next_segregating_site++;
            } else {
                if (n == sample_size) {
                    break;
                }
                for (int i = 0; i < next_segregating_site; i++) {
                    col_sums[i] += gsl_spmatrix_get(tmp, child, i);
                }
                row_counts[child]++;
                n++;
            }
        }
        data = gsl_matrix_alloc(tmp->size1, tmp->size2);
        tmp_data = gsl_matrix_alloc(tmp->size1, tmp->size2);
        gsl_spmatrix_sp2d(data, tmp);
        gsl_spmatrix_free(tmp);
        return;
    }

    void simplify(gsl_matrix *tmp_data, std::vector<int> &tmp_row_counts,
                  std::vector<int> &tmp_col_sums, std::vector<int> &active) {
        bool changed = false;
        int diff = 0, ub_j = 0, lb_i = 0;
        for (unsigned int i = 0; i < tmp_row_counts.size(); i++) {
            changed = false;
            if (tmp_row_counts[i] == 1) {
                for (unsigned int j = 0; j < tmp_col_sums.size(); j++) {
                    if (gsl_matrix_get(tmp_data, i, j) == 1 &&
                        tmp_col_sums[j] == 1) {
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
                            diff += std::abs(gsl_matrix_get(tmp_data, i, k) -
                                             gsl_matrix_get(tmp_data, j, k));
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
                                    active.erase(active.begin() + lb_i + 1,
                                                 active.begin() + lb_i + 2);
                                } else {
                                    active.erase(active.begin() + lb_i,
                                                 active.begin() + lb_i + 1);
                                }
                            }
                            tmp_row_counts[j] += tmp_row_counts[i];
                            tmp_row_counts[i] = 0;
                            for (unsigned int k = 0; k < tmp_col_sums.size();
                                 k++) {
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
        double sim_time = 0;
        gsl_matrix_memcpy(tmp_data, data);
        std::vector<int> tmp_row_counts = row_counts, tmp_col_sums = col_sums;
        std::vector<int> active(sample_size);
        for (int i = 0; i < sample_size; i++) {
            active[i] = i;
        }
        simplify(tmp_data, tmp_row_counts, tmp_col_sums, active);
        int next_parent = sample_size;
        std::vector<double> p(tmp_row_counts.size(), 0);
        double sum_p = 0, ub = 0, coin = 0, time_increment = 0;
        int least_index = 0, c1 = 0, c2 = 0, row = 0;
        while (active.size() > 1) {
            sum_p = 0;
            for (unsigned int i = 0; i < tmp_row_counts.size(); i++) {
                if (tmp_row_counts[i] > 1) {
                    p[i] = gsl_sf_choose(tmp_row_counts[i], 2);
                } else {
                    p[i] = 0;
                }
                sum_p += p[i];
            }
            time_increment =
                gsl_ran_exponential(gen, 1 / gsl_sf_choose(active.size(), 2));
            row = 0;
            ub = p[0];
            coin = gsl_rng_uniform(gen) * sum_p;
            while (coin > ub) {
                row++;
                ub += p[row];
            }
            sim_time += time_increment;
            event_times[next_parent - sample_size] = sim_time;
            event_to_parent[next_parent - sample_size] = next_parent;
            least_index = 0;
            for (int i = 0; i < row; i++) {
                least_index += tmp_row_counts[i];
            }
            c1 = least_index + gsl_rng_uniform_int(gen, tmp_row_counts[row]);
            do {
                c2 =
                    least_index + gsl_rng_uniform_int(gen, tmp_row_counts[row]);
            } while (c1 == c2);
            parent[active[c1]] = next_parent;
            parent[active[c2]] = next_parent;
            child_1[next_parent] = active[c1];
            child_2[next_parent] = active[c2];
            node_times[next_parent] = sim_time;
            tmp_row_counts[row]--;
            active[std::min(c1, c2)] = next_parent;
            active.erase(active.begin() + std::max(c1, c2));
            next_parent++;
            for (unsigned int j = 0; j < tmp_col_sums.size(); j++) {
                tmp_col_sums[j] -= gsl_matrix_get(tmp_data, row, j);
            }
            if (tmp_row_counts[row] == 1) {
                simplify(tmp_data, tmp_row_counts, tmp_col_sums, active);
            }
        }
        for (int i = 0; i < 2 * sample_size - 1; i++) {
            nmut[i] = number_of_mutations(i);
        }
        return;
    }

    int number_of_mutations(const int i) const {
        int ret = 0;
        if (i < sample_size) {
            for (unsigned int j = 0; j < col_sums.size(); j++) {
                if (gsl_matrix_get(data, sample_to_row(i), j) == 1 &&
                    col_sums[j] == 1) {
                    ret++;
                }
            }
        } else {
            unsigned int ind = 0;
            int leaf_sum = 0;
            std::vector<int> leaves(2, 0);
            if (parent[i] != -1) {
                leaves.resize(2);
                leaves[0] = child_1[i];
                leaves[1] = child_2[i];
                ind = 0;
                while (ind < leaves.size()) {
                    if (child_1[leaves[ind]] == -1) {
                        ind++;
                    } else {
                        leaves.push_back(child_2[leaves[ind]]);
                        leaves[ind] = child_1[leaves[ind]];
                    }
                }
                for (unsigned int j = 0; j < col_sums.size(); j++) {
                    leaf_sum = 0;
                    for (unsigned int k = 0; k < leaves.size(); k++) {
                        leaf_sum +=
                            gsl_matrix_get(data, sample_to_row(leaves[k]), j);
                    }
                    if (leaf_sum == (int)leaves.size() &&
                        (int)leaves.size() == col_sums[j]) {
                        ret++;
                    }
                }
            }
        }
        return ret;
    }

    double log_likelihood_times() const {
        double ret =
            -exp(gsl_sf_lnchoose(sample_size, 2) + log(event_times[0]));
        for (int i = 1; i < sample_size - 1; i++) {
            ret -= exp(gsl_sf_lnchoose(sample_size - i, 2) +
                       log(event_times[i] - event_times[i - 1]));
        }
        return ret;
    }

    int sample_to_row(const int sample_index) const {
        int row = 0;
        int total = row_counts[0];
        while (sample_index + 1 > total) {
            row++;
            total += row_counts[row];
        }
        return row;
    }

    double log_likelihood_types() const {
        double ret = 0;
        for (int i = 0; i < 2 * sample_size - 1; i++) {
            if (parent[i] != -1) {
                ret += nmut[i] * (log(mutation_rate) +
                                  log(node_times[parent[i]] - node_times[i])) -
                       mutation_rate * (node_times[parent[i]] - node_times[i]) -
                       gsl_sf_lnfact(nmut[i]);
            }
        }
        return ret;
    }

    void detach_reattach(const int detach, const int attach_above,
                         const double old_parent_time,
                         const double new_parent_time) {
        int mother = parent[detach];
        int sibling = child_2[mother];
        if (detach == sibling) {
            sibling = child_1[mother];
        }
        node_times[mother] = new_parent_time;
        if (attach_above != mother && attach_above != sibling) {
            parent[sibling] = parent[mother];
            if (sibling == child_1[mother]) {
                child_1[mother] = attach_above;
            } else {
                child_2[mother] = attach_above;
            }
            if (parent[mother] != -1) {
                if (child_1[parent[mother]] == mother) {
                    child_1[parent[mother]] = sibling;
                } else {
                    child_2[parent[mother]] = sibling;
                }
            }
            parent[mother] = parent[attach_above];
            if (parent[attach_above] != -1) {
                if (child_1[parent[attach_above]] == attach_above) {
                    child_1[parent[attach_above]] = mother;
                } else {
                    child_2[parent[attach_above]] = mother;
                }
            }
            parent[attach_above] = mother;
        }
        nmut[sibling] = number_of_mutations(sibling);
        nmut[attach_above] = number_of_mutations(attach_above);
        nmut[mother] = number_of_mutations(mother);
        event_to_parent.erase(
            event_to_parent.begin() +
            std::distance(event_times.begin(),
                          std::lower_bound(event_times.begin(),
                                           event_times.end(),
                                           old_parent_time)));
        event_times.erase(std::lower_bound(event_times.begin(),
                                           event_times.end(), old_parent_time));
        event_to_parent.insert(
            event_to_parent.begin() +
                std::distance(event_times.begin(),
                              std::upper_bound(event_times.begin(),
                                               event_times.end(),
                                               new_parent_time)),
            mother);
        event_times.insert(std::upper_bound(event_times.begin(),
                                            event_times.end(), new_parent_time),
                           new_parent_time);
        return;
    }

    int mutations_in_between(int i, int j) const {
        int ret = 0;
        while (i != j) {
            if (j == -1 || node_times[i] <= node_times[j]) {
                ret += nmut[i];
                i = parent[i];
            } else {
                ret += nmut[j];
                j = parent[j];
            }
        }
        return ret;
    }

    int metropolis_topology() {
        int ret = 1;
        int mother = sample_size + gsl_rng_uniform_int(gen, sample_size - 1);
        int detach, sibling;
        if (gsl_rng_uniform(gen) < 0.5) {
            sibling = child_2[mother];
            detach = child_1[mother];
        } else {
            sibling = child_1[mother];
            detach = child_2[mother];
        }
        int attach_above = gsl_rng_uniform_int(gen, 2 * sample_size - 1);
        while (attach_above == detach) {
            attach_above = gsl_rng_uniform_int(gen, 2 * sample_size - 1);
        }
        double log_alpha = 0;
        if (parent[attach_above] != -1 && child_1[detach] != -1 &&
            !(node_times[parent[attach_above]] > node_times[detach])) {
            ret = 0;
        } else if (mutations_in_between(mother, parent[attach_above]) > 0) {
            ret = 0;
        } else {
            double old_parent_time = node_times[mother];
            log_alpha = -ll_times - ll_types;
            double new_parent_time =
                fmax(node_times[attach_above], node_times[detach]);
            double inc = 0;
            if (parent[attach_above] == -1) {
                inc = gsl_ran_exponential(gen, 1);
                log_alpha += inc;
            } else {
                inc = gsl_rng_uniform(gen) *
                      (node_times[parent[attach_above]] - new_parent_time);
                log_alpha +=
                    log(node_times[parent[attach_above]] - new_parent_time);
            }
            new_parent_time += inc;
            if (parent[mother] == -1) {
                log_alpha -= node_times[mother] -
                             fmax(node_times[sibling], node_times[detach]);
            } else {
                log_alpha -= log(node_times[parent[mother]] -
                                 fmax(node_times[sibling], node_times[detach]));
            }
            double old_ll_times = ll_times;
            double old_ll_types = ll_types;
            detach_reattach(detach, attach_above, old_parent_time,
                            new_parent_time);
            ll_times = log_likelihood_times();
            ll_types = log_likelihood_types();
            log_alpha += ll_times + ll_types;
            if (log(gsl_rng_uniform_pos(gen)) > log_alpha) {
                ret = 0;
                ll_times = old_ll_times;
                ll_types = old_ll_types;
                detach_reattach(detach, sibling, new_parent_time,
                                old_parent_time);
            }
        }
        return ret;
    }

    int metropolis_times(const double sd_numerator) {
        int ret = 1;
        double log_alpha = -ll_times - ll_types;
        std::vector<double> old_event_times = event_times;
        std::vector<double> old_node_times = node_times;
        std::vector<int> old_event_to_parent = event_to_parent;
        double lb, mu, seed, sd;
        for (int i = 0; i < sample_size - 1; i++) {
            sd = sd_numerator / sqrt((sample_size - 1) * (sample_size - i) *
                                     (sample_size - i - 1));
            lb = fmax(node_times[child_1[old_event_to_parent[i]]],
                      node_times[child_2[old_event_to_parent[i]]]);
            mu = old_event_times[i];
            seed = gsl_rng_uniform(gen);
            seed += (1 - seed) * gsl_cdf_gaussian_P(lb - mu, sd);
            seed = mu + gsl_cdf_gaussian_Pinv(seed, sd);
            detach_reattach(child_1[old_event_to_parent[i]],
                            child_2[old_event_to_parent[i]], mu, seed);
            log_alpha -= log(gsl_ran_gaussian_pdf(seed - mu, sd)) -
                         log(1 - gsl_cdf_gaussian_P(lb - mu, sd));
        }
        for (int i = 0; i < sample_size - 1; i++) {
            sd = sd_numerator / sqrt((sample_size - 1) * (sample_size - i) *
                                     (sample_size - i - 1));
            lb = fmax(old_node_times[child_1[event_to_parent[i]]],
                      old_node_times[child_2[event_to_parent[i]]]);
            mu = event_times[i];
            log_alpha += log(gsl_ran_gaussian_pdf(
                             old_node_times[event_to_parent[i]] - mu, sd)) -
                         log(1 - gsl_cdf_gaussian_P(lb - mu, sd));
        }
        double old_ll_times = ll_times;
        double old_ll_types = ll_types;
        ll_times = log_likelihood_times();
        ll_types = log_likelihood_types();
        log_alpha += ll_times + ll_types;
        if (log(gsl_rng_uniform_pos(gen)) > log_alpha) {
            ret = 0;
            ll_times = old_ll_times;
            ll_types = old_ll_types;
            event_times = old_event_times;
            node_times = old_node_times;
            event_to_parent = old_event_to_parent;
        }
        return ret;
    }

    int metropolis_mutation_rate(const double sd) {
        int ret = 1;
        double old_mutation_rate = mutation_rate;
        double log_alpha = -ll_types;
        // sd divided by two because the algorithm works with mutation_rate / 2
        mutation_rate += gsl_ran_gaussian_ziggurat(gen, sd / 2);
        if (mutation_rate < 0.0) {
            mutation_rate = -mutation_rate;
        }
        double old_ll_types = ll_types;
        ll_types = log_likelihood_types();
        log_alpha += ll_types;
        if (log(gsl_rng_uniform(gen)) > log_alpha) {
            ret = 0;
            mutation_rate = old_mutation_rate;
            ll_types = old_ll_types;
        }
        return ret;
    }

    double tree_height() const { return event_times.back(); }

    double mutation_rate, ll_times, ll_types;
    int sample_size;
    std::vector<int> row_counts, col_sums, parent, child_1, child_2,
        event_to_parent, nmut;
    std::vector<double> node_times, event_times;
    gsl_matrix *data;
    gsl_matrix *tmp_data;
    gsl_rng *gen;
};

#endif
