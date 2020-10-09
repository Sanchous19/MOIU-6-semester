#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

using vector = std::vector<double>;


const double INF = 1e9;
const double EPS = 1e-12;


template <typename T>
std::istream& operator>>(std::istream& in, std::vector<T>& v)
{
	for (T& num : v)
		in >> num;
	return in;
}


template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v)
{
	for (const T& num : v)
		out << num << " ";
	return out;
}


vector operator*(const double& value, const vector& v)
{
	vector result(v.size());
	for (size_t i = 0; i < v.size(); i++)
		result[i] = value * v[i];
	return result;
}


void operator-=(vector& v1, const vector& v2)
{
	if (v1.size() != v2.size())
		throw std::invalid_argument("Can't difference two vectors, there sizes are not equal");

	for (size_t i = 0; i < v1.size(); i++)
		v1[i] -= v2[i];
}


void operator/=(vector& v, const double& value)
{
	for (double& num : v)
		num /= value;
}


class Matrix
{
public:
	Matrix(const size_t& row_num, const size_t& col_num)
		: row_num_(row_num), col_num_(col_num), data_(col_num, vector(row_num)) {}

	void MakeIdentity()
	{
		if (row_num_ != col_num_)
			throw std::invalid_argument("Can't make identity matrix, the number of rows is not equal to the number of column");

		for (size_t i = 0; i < row_num_; i++)
			for (size_t j = 0; j < col_num_; j++)
				data_[j][i] = (i == j) ? 1 : 0;
	}

	Matrix Inverse() const
	{
		if (row_num_ != col_num_)
			throw std::invalid_argument("Can't get inverse matrix, the number of rows is not equal to the number of column");

		size_t n = row_num_;
		Matrix base = *this;
		Matrix inv(n, n);
		inv.MakeIdentity();

		for (size_t i = 0; i < n; i++)
		{
			size_t col = i;
			for (size_t j = i + 1; j < n; j++)
			{
				if (abs(base[j][i]) > abs(base[col][i]))
					col = j;
			}
			swap(base[col], base[i]);
			swap(inv[col], inv[i]);

			double diag_elem = base[i][i];
			base[i] /= diag_elem;
			inv[i] /= diag_elem;

			for (size_t j = i + 1; j < n; j++)
			{
				double coeff = base[j][i];
				base[j] -= coeff * base[i];
				inv[j] -= coeff * inv[i];
			}
		}
		for (size_t i = 0; i < n; i++)
		{
			for (size_t j = 0; j < n - 1 - i; j++)
			{
				double coeff = base[j][n - 1 - i];
				base[j] -= coeff * base[n - 1 - i];
				inv[j] -= coeff * inv[n - 1 - i];
			}
		}
		return inv;
	}

	void Extend(Matrix& m)
	{
		data_.reserve(m.col_num_);
		for (vector& col : m)
			data_.push_back(std::move(col));
		col_num_ += m.col_num_;
	}

	void Resize(const size_t& size)
	{
		data_.resize(size);
	}

	void DeleteRow(const size_t& n)
	{
		for (vector& row : data_)
			row.erase(row.begin() + n);
		row_num_--;
	}

	size_t RowSize() const
	{
		return row_num_;
	}

	size_t ColumnSize() const
	{
		return col_num_;
	}

	void SetRowSize(const size_t& value)
	{
		row_num_ = value;
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
	friend Matrix operator*(const Matrix& m1, const Matrix& m2);
	friend vector operator*(const vector& v, const Matrix& m);
	friend vector operator*(const Matrix& m, const vector& v);

private:
	std::vector<vector> data_;
	size_t row_num_, col_num_;
};


std::istream& operator>>(std::istream& in, Matrix& m)
{
	for (size_t i = 0; i < m.row_num_; i++)
		for (size_t j = 0; j < m.col_num_; j++)
			in >> m[j][i];
	return in;
}


std::ostream& operator<<(std::ostream& out, const Matrix& m)
{
	for (size_t i = 0; i < m.row_num_; i++)
	{
		for (size_t j = 0; j < m.col_num_; j++)
			out << m[j][i] << " ";
		if (i + 1 != m.row_num_)
			out << "\n";
	}
	return out;
}


Matrix operator*(const Matrix& m1, const Matrix& m2)
{
	if (m1.col_num_ != m2.row_num_)
		throw std::invalid_argument("Can't multiply two matrices, the number of columns of the first matrix "
			"is not equal to the number of rows of the second matrix");

	size_t row_num = m1.row_num_, col_num = m2.col_num_, sum_num = m1.col_num_;
	Matrix result(row_num, col_num);
	for (size_t i = 0; i < row_num; i++)
		for (size_t j = 0; j < col_num; j++)
		{
			result[j][i] = 0;
			for (size_t k = 0; k < sum_num; k++)
				result[j][i] += m1[k][i] * m2[j][k];
		}

	return result;
}


vector operator*(const Matrix& m, const vector& v)
{
	if (m.col_num_ != v.size())
		throw std::invalid_argument("Can't multiply matrix and vector, the number of columns "
			"of the matrix is not equal to the size of vector");

	vector result(m.row_num_, 0);
	for (size_t i = 0; i < m.row_num_; i++)
		for (size_t j = 0; j < m.col_num_; j++)
			result[i] += m[j][i] * v[j];

	return result;
}


vector operator*(const vector& v, const Matrix& m)
{
	if (v.size() != m.row_num_)
		throw std::invalid_argument("Can't multiply vector and matrix, the size of vector "
			"is not equal the number of rows of matrix");

	vector result(m.col_num_, 0);
	for (size_t i = 0; i < m.col_num_; i++)
		for (size_t j = 0; j < m.row_num_; j++)
			result[i] += v[j] * m[i][j];

	return result;
}


double operator*(const vector& v1, const vector& v2)
{
	if (v1.size() != v2.size())
		throw std::invalid_argument("Can't multiply two vectors, there sizes are not equal");

	double result = 0;
	for (size_t i = 0; i < v1.size(); i++)
		result += v1[i] * v2[i];
	return result;
}


void TransformInverseMatrix(Matrix& B, const vector& zv, const size_t& s)
{
	size_t n = B.RowSize();
	vector sv(n);
	for (size_t i = 0; i < n; i++)
		sv[i] = B[i][s];

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			if (i != s)
				B[j][i] -= sv[j] * zv[i] / zv[s];
			else
				B[j][i] = sv[j] / zv[s];
		}
	}
}


bool SimplexAlgorithm(const size_t& n, const size_t& m, vector& c, Matrix& A, vector& x, std::vector<size_t>& jb, Matrix* pB = nullptr)
{
	bool is_has_solve = true;

	std::set<size_t> jn;
	for (size_t i = 0; i < n; i++)
		jn.insert(i);

	Matrix Ab(0, m);
	vector cb(m);
	for (size_t i = 0; i < m; i++)
	{
		size_t j = jb[i];
		Ab[i] = std::move(A[j]);
		cb[i] = c[j];
		jn.erase(j);
	}
	Ab.SetRowSize(m);

	Matrix B = Ab.Inverse();

	while (true)
	{
		vector u = cb * B;
		double min_delta = 0;
		size_t j0 = 0;
		for (const size_t& j : jn)
		{
			double delta = u * A[j] - c[j];
			if (delta < min_delta)
			{
				min_delta = delta;
				j0 = j;
			}
		}
		if (min_delta >= -EPS)
			break;

		vector zv = B * A[j0];
		size_t s = 0;
		double min_teta = INF;
		for (size_t i = 0; i < m; i++)
		{
			size_t j = jb[i];
			double z = zv[i];
			double teta = z > EPS ? x[j] / z : INF;
			if (teta < min_teta)
			{
				min_teta = teta;
				s = i;
			}
		}
		if (min_teta == INF)
		{
			is_has_solve = false;
			break;
		}

		size_t js = jb[s];
		A[js] = std::move(Ab[s]);
		Ab[s] = std::move(A[j0]);
		cb[s] = c[j0];
		jb[s] = j0;
		jn.erase(j0);
		jn.insert(js);

		for (size_t i = 0; i < m; i++)
		{
			size_t j = jb[i];
			x[j] -= min_teta * zv[i];
		}
		x[j0] = min_teta;
		x[js] = 0;

		TransformInverseMatrix(B, zv, s);
	}

	for (size_t i = 0; i < m; i++)
	{
		size_t j = jb[i];
		A[j] = std::move(Ab[i]);
	}

	if (pB)
		*pB = std::move(B);

	return is_has_solve;
}


Matrix GetNewInverseMatrix(Matrix& A, const std::vector<size_t>& jb)
{
	size_t n = jb.size();
	Matrix Ab(0, n);

	for (size_t i = 0; i < n; i++)
	{
		size_t j = jb[i];
		Ab[i] = std::move(A[j]);
	}
	Ab.SetRowSize(n);

	Matrix B = Ab.Inverse();
	for (size_t i = 0; i < n; i++)
	{
		size_t j = jb[i];
		A[j] = std::move(Ab[i]);
	}

	return B;
}


bool BuildBaseSolution(const size_t& n, size_t& m, Matrix& A, vector& b, vector& x, std::vector<size_t>& jb)
{
	x.assign(n + m, 0);
	jb.resize(m);
	vector c(n + m, 0);
	for (size_t i = n; i < n + m; i++)
	{
		x[i] = b[i - n];
		jb[i - n] = i;
		c[i] = -1;
	}

	Matrix E(m, m);
	E.MakeIdentity();
	A.Extend(E);

	Matrix* pB = new Matrix(0, 0);
	SimplexAlgorithm(n, m, c, A, x, jb, pB);

	for (size_t i = n; i < n + m; i++)
	{
		if (x[i])
		{
			A.Resize(m);
			return false;
		}
	}

	std::set<size_t> jn;
	for (size_t i = 0; i < n; i++)
		jn.insert(i);

	std::vector<size_t> ub;
	for (size_t k = 0; k < m; k++)
	{
		size_t j = jb[k];
		if (j >= n)
			ub.push_back(k);
		else
			jn.erase(j);
	}
	std::sort(ub.begin(), ub.end());

	Matrix& B = *pB;
	std::set<size_t> delete_rows;
	for (size_t& k : ub)
	{
		k -= delete_rows.size();
		size_t jk = jb[k], j0 = 0;
		double max_alpha = 0;
		for (const size_t& j : jn)
		{
			double alpha = 0;
			for (size_t i = 0; i < jb.size(); i++)
				alpha += B[i][k] * A[j][i];

			if (abs(alpha) >= abs(max_alpha))
			{
				max_alpha = alpha;
				j0 = j;
			}
		}

		if (max_alpha)
		{
			vector zv = B * A[j0];
			TransformInverseMatrix(B, zv, k);
			jb[k] = j0;
			jn.erase(j0);
		}
		else
		{
			size_t row = jk - n - std::distance(delete_rows.begin(), delete_rows.lower_bound(jk));
			jb.erase(jb.begin() + k);
			A.DeleteRow(row);
			B = GetNewInverseMatrix(A, jb);
			delete_rows.insert(jk);
		}
	}
	m -= delete_rows.size();

	A.Resize(n);
	x.resize(n);
	return true;
}


void PrepareMatrix(const size_t& n, const size_t& m, Matrix& A, vector& b)
{
	for (size_t i = 0; i < m; i++)
	{
		if (b[i] < 0)
		{
			b[i] = -b[i];
			for (size_t j = 0; j < n; j++)
				A[j][i] = -A[j][i];
		}
	}
}


int main()
{
	std::ios_base::sync_with_stdio(false);
	std::cin.tie(0);

	size_t n, m;
	std::cin >> m >> n;

	Matrix A(m, n);
	std::cin >> A;

	vector b(m);
	std::cin >> b;

	vector c(n);
	std::cin >> c;

	try
	{
		PrepareMatrix(n, m, A, b);

		vector x;
		std::vector<size_t> jb;
		bool is_has_solve = BuildBaseSolution(n, m, A, b, x, jb);

		if (!is_has_solve)
		{
			std::cout << "No solution";
		}
		else
		{
			bool is_bounded = SimplexAlgorithm(n, m, c, A, x, jb);

			if (is_bounded)
			{
				std::cout << "Bounded\n";
				std::cout << std::fixed << std::setprecision(9);
				std::cout << x << "\n";
			}
			else
			{
				std::cout << "Unbounded";
			}
		}
	}
	catch (std::exception & e)
	{
		std::cout << e.what();
	}

	return 0;
}
