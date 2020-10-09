#include <algorithm>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <vector>

using vector = std::vector<int>;
using cell = std::pair<size_t, size_t>;
using cell_set = std::set<cell>;
using cell_vector = std::vector<cell>;


template <typename T>
std::istream& operator>>(std::istream& in, std::vector<T>& v)
{
	for (T& num : v)
		in >> num;
	return in;
}


class Matrix
{
public:
	Matrix(const size_t& row_num, const size_t& col_num)
		: row_num_(row_num), col_num_(col_num), data_(row_num, std::vector<int>(col_num)) {}

	void MakeZero()
	{
		for (vector& row : data_)
			for (int& val : row)
				val = 0;
	}

	void PushRow()
	{
		data_.push_back(vector(col_num_, 0));
		row_num_++;
	}

	void PushColumn()
	{
		for (vector& row : data_)
			row.push_back(0);
		col_num_++;
	}

	void PopRow()
	{
		data_.pop_back();
		row_num_--;
	}

	void PopColumn()
	{
		for (vector& row : data_)
			row.pop_back();
		col_num_--;
	}

	size_t RowSize() const
	{
		return row_num_;
	}

	size_t ColumnSize() const
	{
		return col_num_;
	}

	typename std::vector<vector>::iterator begin()
	{
		return data_.begin();
	}

	typename std::vector<vector>::const_iterator begin() const
	{
		return data_.begin();
	}

	typename std::vector<vector>::iterator end()
	{
		return data_.end();
	}

	typename std::vector<vector>::const_iterator end() const
	{
		return data_.end();
	}

	vector& operator[](const size_t& index)
	{
		return data_[index];
	}

	const vector& operator[](const size_t& index) const
	{
		return data_[index];
	}

	friend std::istream& operator>>(std::istream& in, Matrix& m);
	friend std::ostream& operator<<(std::ostream& out, const Matrix& m);

private:
	std::vector<vector> data_;
	size_t row_num_, col_num_;
};


std::istream& operator>>(std::istream& in, Matrix& m)
{
	for (vector& row : m)
		for (int& val : row)
			in >> val;
	return in;
}


std::ostream& operator<<(std::ostream& out, const Matrix& m)
{
	for (size_t i = 0; i < m.row_num_; i++)
	{
		for (const int& val : m[i])
			out << val << " ";
		if (i + 1 != m.row_num_)
			out << "\n";
	}
	return out;
}


int Sum(const vector& v)
{
	int res = 0;
	for (const int& item : v)
		res += item;
	return res;
}


cell_set ConstructBase(vector a, vector b, Matrix& X)
{
	X.MakeZero();
	cell_set base;
	for (size_t i = 0, j = 0, k = 0; k < (a.size() + b.size() - 1); k++)
	{
		int min_val = std::min(a[i], b[j]);
		a[i] -= min_val;
		b[j] -= min_val;
		X[i][j] = min_val;
		base.insert({ i, j });

		if (j < b.size() - 1 && !b[j])
			j++;
		else if (i < a.size() - 1 && !a[i])
			i++;
	}
	return base;
}


bool GetNewBaseCell(cell& new_b, const std::vector<cell_set>& row_base, const std::vector<cell_set>& column_base, const Matrix& C)
{
	size_t m = C.RowSize(), n = C.ColumnSize();

	vector u(m, 0), v(n, 0);
	std::vector<bool> ub(m, false), vb(n, false);

	std::queue<size_t> q;
	q.push(0);
	ub[0] = true;
	while (!q.empty())
	{
		size_t line = q.front();
		q.pop();
		if (line < m)
		{
			for (const cell& b : row_base[line])
			{
				size_t column = b.second;
				if (!vb[column])
				{
					v[column] = C[line][column] - u[line];
					q.push(column + m);
					vb[column] = true;
				}
			}
		}
		else
		{
			line -= m;
			for (const cell& b : column_base[line])
			{
				size_t row = b.first;
				if (!ub[row])
				{
					u[row] = C[row][line] - v[line];
					q.push(row);
					ub[row] = true;
				}
			}
		}
	}

	int min_value = 0;
	for (size_t i = 0; i < m; i++)
		for (size_t j = 0; j < n; j++)
		{
			int value = C[i][j] - u[i] - v[j];
			if (value < min_value)
			{
				min_value = value;
				new_b = { i, j };
			}
		}

	return min_value < 0;
}


void ChangeBase(cell_set& base, std::vector<cell_set>& row_base, std::vector<cell_set>& column_base, vector& row_base_cnt, vector& column_base_cnt, cell& new_b, Matrix& X)
{
	base.insert(new_b);
	row_base[new_b.first].insert(new_b);
	column_base[new_b.second].insert(new_b);
	row_base_cnt[new_b.first]++;
	column_base_cnt[new_b.second]++;

	vector copy_row_base_cnt = row_base_cnt, copy_column_base_cnt = column_base_cnt;
	std::vector<bool> is_delete(base.size(), false);
	bool is_circle = false;
	while (!is_circle)
	{
		is_circle = true;
		size_t i = 0;
		for (const cell& b : base)
		{
			i++;
			if (is_delete[i - 1])
				continue;

			size_t row = b.first, column = b.second;
			if (copy_row_base_cnt[row] == 1 || copy_column_base_cnt[column] == 1)
			{
				copy_row_base_cnt[row]--;
				copy_column_base_cnt[column]--;
				is_delete[i - 1] = true;
				is_circle = false;
			}
		}
	}

	cell_vector graph_nodes;
	size_t i = 0;
	for (const cell& b : base)
	{
		i++;
		if (!is_delete[i - 1])
			graph_nodes.push_back(b);
	}

	std::map<cell, cell_vector> graph;
	std::sort(graph_nodes.begin(), graph_nodes.end());
	for (size_t i = 1; i < graph_nodes.size(); i++)
		if (graph_nodes[i].first == graph_nodes[i - 1].first)
		{
			graph[graph_nodes[i]].push_back(graph_nodes[i - 1]);
			graph[graph_nodes[i - 1]].push_back(graph_nodes[i]);
		}

	std::sort(graph_nodes.begin(), graph_nodes.end(),
		[](const cell& node1, const cell& node2) {
		return node1.second < node2.second;
	});
	for (size_t i = 1; i < graph_nodes.size(); i++)
		if (graph_nodes[i].second == graph_nodes[i - 1].second)
		{
			graph[graph_nodes[i]].push_back(graph_nodes[i - 1]);
			graph[graph_nodes[i - 1]].push_back(graph_nodes[i]);
		}

	cell_vector way;
	int min_x = 1e9;
	bool is_positive = true;
	cell prev_node = new_b, next_node = new_b;
	do
	{
		if (!is_positive)
			min_x = std::min(min_x, X[next_node.first][next_node.second]);
		is_positive = !is_positive;
		way.push_back(next_node);
		for (const cell& node : graph[next_node])
		{
			if (node != prev_node)
			{
				prev_node = next_node;
				next_node = node;
				break;
			}
		}
	} while (next_node != new_b);

	cell delete_b;
	int min_new_x = 1e9;
	is_positive = true;
	for (const cell& node : way)
	{
		if (is_positive)
			X[node.first][node.second] += min_x;
		else
		{
			X[node.first][node.second] -= min_x;
			if (X[node.first][node.second] < min_new_x)
			{
				delete_b = node;
				min_new_x = X[node.first][node.second];
			}
		}
		is_positive = !is_positive;
	}

	base.erase(delete_b);
	row_base[delete_b.first].erase(delete_b);
	column_base[delete_b.second].erase(delete_b);
	row_base_cnt[delete_b.first]--;
	column_base_cnt[delete_b.second]--;
}


Matrix Solve(const vector& a, const vector& b, const Matrix& C)
{
	size_t m = a.size(), n = b.size();

	Matrix X(m, n);
	cell_set base = ConstructBase(a, b, X);
	std::vector<cell_set> row_base(m), column_base(n);
	vector row_base_cnt(m), column_base_cnt(n);
	for (const cell& b : base)
	{
		row_base[b.first].insert(b);
		column_base[b.second].insert(b);
		row_base_cnt[b.first]++;
		column_base_cnt[b.second]++;
	}

	while (true)
	{
		cell new_b;
		if (!GetNewBaseCell(new_b, row_base, column_base, C))
			break;
		ChangeBase(base, row_base, column_base, row_base_cnt, column_base_cnt, new_b, X);
	}

	return X;
}


int main()
{
	std::ios_base::sync_with_stdio(false);
	std::cin.tie(0);

	size_t n, m;
	std::cin >> m >> n;

	Matrix C(m, n);
	std::cin >> C;

	vector a(m);
	std::cin >> a;

	vector b(n);
	std::cin >> b;

	int sum_a = Sum(a), sum_b = Sum(b);
	if (sum_a > sum_b)
	{
		C.PushColumn();
		b.push_back(sum_a - sum_b);
	}
	else if (sum_a < sum_b)
	{
		C.PushRow();
		a.push_back(sum_b - sum_a);
	}

	Matrix X = Solve(a, b, C);

	if (X.ColumnSize() != n)
		X.PopColumn();
	else if (X.RowSize() != m)
		X.PopRow();
	std::cout << X;

	return 0;
}
