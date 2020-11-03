#ifndef SIM
#define SIM

#include <algorithm>
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
        : mutation_rate(), ll_times(), ll_types(), sample_size(0), nsites(),
          row_counts(), parent(), child_1(), child_2(), event_to_parent(),
          node_times(), event_times(), L() {
        gen = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(gen, time(NULL) * getpid());
        std::ifstream file;
        file.open(filename);
        std::string line, token;
        gsl_spmatrix *tmp = gsl_spmatrix_alloc(1, 1);
        int l = 0;
        std::vector<int> tmp_vec;
        while (getline(file, line)) {
            tmp_vec.clear();
            std::stringstream iss;
            iss << line;
            while (getline(iss, token, ' ')) {
                tmp_vec.push_back(atoi(token.c_str()));
            }
            for (unsigned int i = 0; i < tmp_vec.size() - 1; i++) {
                if (tmp_vec[i] == 1) {
                    gsl_spmatrix_set(tmp, l, i, 1);
                }
            }
            row_counts.push_back(tmp_vec.back());
            sample_size += tmp_vec.back();
            l++;
        }
        nsites = tmp_vec.size() - 1;
        parent.resize(2 * sample_size - 1, -1);
        child_1.resize(2 * sample_size - 1, -1);
        child_2.resize(2 * sample_size - 1, -1);
        event_to_parent.resize(sample_size - 1, -1);
        node_times.resize(2 * sample_size - 1, 0.0);
        event_times.resize(sample_size - 1, 0.0);
        std::vector<double> L_tmp(2);
        std::vector<std::vector<double>> L_full(2 * sample_size - 1, L_tmp);
        L.resize(nsites, L_full);
        data = gsl_matrix_alloc(tmp->size1, tmp->size2);
        gsl_spmatrix_sp2d(data, tmp);
        gsl_spmatrix_free(tmp);
        mutation_rate = watterson_estimator() / 2;
        sample_tree();
        ll_times = log_likelihood_times();
        ll_types = log_likelihood_types();
    }

    Sim(const int n, const int m, const double theta)
        : mutation_rate(theta / 2), ll_times(), ll_types(), sample_size(n),
          nsites(m), row_counts(), parent(2 * n - 1, -1),
          child_1(2 * n - 1, -1), child_2(2 * n - 1, -1),
          event_to_parent(n - 1, -1), node_times(n - 1, 0.0),
          event_times(n - 1, 0.0), L() {
        gen = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(gen, time(NULL) * getpid());
    }

    ~Sim() {
        gsl_rng_free(gen);
        gsl_matrix_free(data);
    }

    double watterson_estimator() const {
        double harmonic = 1;
        for (int i = 2; i < sample_size; i++) {
            harmonic += 1.0 / i;
        }
        int seg = 0;
        int anc;
        for (int i = 0; i < nsites; i++) {
            anc = gsl_matrix_get(data, 0, i);
            for (unsigned int j = 1; j < row_counts.size(); j++) {
                if (gsl_matrix_get(data, j, i) != anc) {
                    seg++;
                    break;
                }
            }
        }
        return seg / harmonic;
    }

    void generate_data() {
        gsl_spmatrix *tmp = gsl_spmatrix_alloc(1, nsites);
        for (int i = 0; i < nsites; i++) {
            gsl_spmatrix_set(tmp, 0, i, gsl_rng_uniform_int(gen, 2));
        }
        double theta = 2 * mutation_rate;
        int n = 2;
        row_counts.push_back(2);
        int child = 0;
        int ub = row_counts[0];
        double coin = 0;
        while (n < sample_size + 1) {
            child = 0;
            ub = row_counts[0];
            coin = gsl_rng_uniform(gen);
            while (coin > double(ub) / double(n)) {
                child++;
                ub += row_counts[child];
            }
            if (gsl_rng_uniform(gen) < theta / (n - 1 + theta)) {
                if (row_counts[child] > 1) {
                    row_counts[child]--;
                    row_counts.push_back(1);
                    for (int i = 0; i < nsites; i++) {
                        // set to 1 first to make GSL resize the matrix properly
                        // even if the new entry is 0
                        gsl_spmatrix_set(tmp, (int)(row_counts.size() - 1), i,
                                         1);
                        gsl_spmatrix_set(tmp, (int)(row_counts.size() - 1), i,
                                         gsl_spmatrix_get(tmp, child, i));
                    }
                    child = row_counts.size() - 1;
                }
                int site = gsl_rng_uniform_int(gen, nsites);
                gsl_spmatrix_set(tmp, child, site,
                                 ((int)gsl_spmatrix_get(tmp, child, site) + 1) %
                                     2);
                int check;
                for (unsigned int i = 0; i < row_counts.size(); i++) {
                    check = 0;
                    if ((int)i != child) {
                        check = 1;
                        for (int j = 0; j < nsites; j++) {
                            if (gsl_spmatrix_get(tmp, i, j) !=
                                gsl_spmatrix_get(tmp, child, j)) {
                                check = 0;
                                break;
                            }
                        }
                    }
                    if (check == 1) {
                        row_counts[i]++;
                        row_counts[child]--;
                        break;
                    }
                }
            } else {
                if (n == sample_size) {
                    break;
                }
                row_counts[child]++;
                n++;
            }
        }
        int r = 0;
        for (unsigned int i = 0; i < row_counts.size(); i++) {
            if (row_counts[i] > 0) {
                r++;
            }
        }
        data = gsl_matrix_alloc(r, nsites);
        r--;
        for (int i = (int)row_counts.size() - 1; i > -1; i--) {
            if (row_counts[i] > 0) {
                for (int j = 0; j < nsites; j++) {
                    gsl_matrix_set(data, r, j, gsl_spmatrix_get(tmp, i, j));
                }
                r--;
            } else {
                row_counts.erase(row_counts.begin() + i);
            }
        }
        gsl_spmatrix_free(tmp);
        return;
    }

    void sample_tree() {
        double sim_time = 0;
        std::vector<int> active(sample_size);
        for (int i = 0; i < sample_size; i++) {
            active[i] = i;
        }
        int next_parent = sample_size;
        for (int i = 0; i < (int)row_counts.size(); i++) {
            for (int j = 0; j < row_counts[i] - 1; j++) {
                sim_time += gsl_ran_exponential(
                    gen, 1 / gsl_sf_choose(active.size(), 2));
                event_times[next_parent - sample_size] = sim_time;
                event_to_parent[next_parent - sample_size] = next_parent;
                parent[active[i]] = next_parent;
                parent[active[i + 1]] = next_parent;
                child_1[next_parent] = active[i];
                child_2[next_parent] = active[i + 1];
                node_times[next_parent] = sim_time;
                active[i] = next_parent;
                active.erase(active.begin() + i + 1);
                next_parent++;
            }
        }
        for (int i = 0; i < (int)row_counts.size() - 1; i++) {
            sim_time +=
                gsl_ran_exponential(gen, 1 / gsl_sf_choose(active.size(), 2));
            event_times[next_parent - sample_size] = sim_time;
            event_to_parent[next_parent - sample_size] = next_parent;
            parent[active[0]] = next_parent;
            parent[active[1]] = next_parent;
            child_1[next_parent] = active[0];
            child_2[next_parent] = active[1];
            node_times[next_parent] = sim_time;
            active[0] = next_parent;
            active.erase(active.begin() + 1);
            next_parent++;
        }
        for (int i = 0; i < nsites; i++) {
            for (int j = 0; j < sample_size; j++) {
                if (gsl_matrix_get(data, sample_to_row(j), i) == 1) {
                    L[i][j][1] = 0;
                    L[i][j][0] = std::numeric_limits<double>::lowest();
                } else {
                    L[i][j][1] = std::numeric_limits<double>::lowest();
                    L[i][j][0] = 0;
                }
            }
        }
        return;
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

    double log_mutation_transition(const double t, const int i,
                                   const int j) const {
        double ret = (1 - exp(-2 * mutation_rate * t / nsites)) / 2;
        if (i == j) {
            ret += exp(-2 * mutation_rate * t / nsites);
        }
        return log(ret);
    }

    double log_likelihood_types() {
        double ret = 0;
        int p, c;
        double t, a, b;
        for (int i = 0; i < nsites; i++) {
            for (int j = 0; j < sample_size - 1; j++) {
                p = event_to_parent[j];
                c = child_1[p];
                t = node_times[p] - node_times[c];
                a = log_mutation_transition(t, 0, 0) + L[i][c][0];
                b = log_mutation_transition(t, 0, 1) + L[i][c][1];
                L[i][p][0] = fmax(a, b) + log(1 + exp(fmin(a, b) - fmax(a, b)));
                a = log_mutation_transition(t, 1, 0) + L[i][c][0];
                b = log_mutation_transition(t, 1, 1) + L[i][c][1];
                L[i][p][1] = fmax(a, b) + log(1 + exp(fmin(a, b) - fmax(a, b)));
                c = child_2[p];
                t = node_times[p] - node_times[c];
                a = log_mutation_transition(t, 0, 0) + L[i][c][0];
                b = log_mutation_transition(t, 0, 1) + L[i][c][1];
                L[i][p][0] +=
                    fmax(a, b) + log(1 + exp(fmin(a, b) - fmax(a, b)));
                a = log_mutation_transition(t, 1, 0) + L[i][c][0];
                b = log_mutation_transition(t, 1, 1) + L[i][c][1];
                L[i][p][1] +=
                    fmax(a, b) + log(1 + exp(fmin(a, b) - fmax(a, b)));
            }
            p = event_to_parent.back();
            a = L[i][p][0];
            b = L[i][p][1];
            ret += fmax(a, b) + log(1 + exp(fmin(a, b) - fmax(a, b))) - log(2);
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

    int metropolis_topology() {
        int ret = 1;
        int detach = gsl_rng_uniform_int(gen, 2 * sample_size - 1);
        while (parent[detach] == -1) {
            detach = gsl_rng_uniform_int(gen, 2 * sample_size - 1);
        }
        int mother = parent[detach];
        int attach_above = gsl_rng_uniform_int(gen, 2 * sample_size - 1);
        while (attach_above == detach) {
            attach_above = gsl_rng_uniform_int(gen, 2 * sample_size - 1);
        }
        double log_alpha = 0;
        int sibling = child_2[mother];
        if (detach == sibling) {
            sibling = child_1[mother];
        }
        if (parent[attach_above] != -1 && child_1[detach] != -1 &&
            !(node_times[parent[attach_above]] > node_times[detach])) {
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
    int sample_size, nsites;
    std::vector<int> row_counts, parent, child_1, child_2, event_to_parent;
    std::vector<double> node_times, event_times;
    std::vector<std::vector<std::vector<double>>> L;
    gsl_matrix *data;
    gsl_rng *gen;
};

#endif
