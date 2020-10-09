#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

using vector = std::vector<double>;
using index_set = std::set<size_t>;


const double INF = std::numeric_limits<double>::infinity();
const double EPS = 1e-6;


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


vector operator-(const vector& v)
{
	vector res(v.size());
	for (size_t i = 0; i < v.size(); i++)
		res[i] = -v[i];
	return res;
}


vector operator*(const double& value, const vector& v)
{
	vector result(v.size());
	for (size_t i = 0; i < v.size(); i++)
		result[i] = value * v[i];
	return result;
}


void operator/=(vector& v, const double& value)
{
	for (double& num : v)
		num /= value;
}


vector operator+(const vector& v1, const vector& v2)
{
	if (v1.size() != v2.size())
		throw std::invalid_argument("Can't difference two vectors, there sizes are not equal");

	size_t n = v1.size();
	vector result(n);
	for (size_t i = 0; i < n; i++)
		result[i] = v1[i] + v2[i];
	return result;
}


void operator+=(vector& v1, const vector& v2)
{
	if (v1.size() != v2.size())
		throw std::invalid_argument("Can't difference two vectors, there sizes are not equal");

	for (size_t i = 0; i < v1.size(); i++)
		v1[i] += v2[i];
}


void operator-=(vector& v1, const vector& v2)
{
	if (v1.size() != v2.size())
		throw std::invalid_argument("Can't difference two vectors, there sizes are not equal");

	for (size_t i = 0; i < v1.size(); i++)
		v1[i] -= v2[i];
}


class Matrix
{
public:
	Matrix(const size_t& row_num, const size_t& col_num)
		: row_num_(row_num), col_num_(col_num), data_(col_num, std::vector<double>(row_num)) {}

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
				if (fabs(base[j][i]) > fabs(base[col][i]))
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

	void PushBackColumn(const vector& v)
	{
		data_.push_back(v);
		col_num_++;
	}

	vector getRow(size_t i)
	{
		vector row(col_num_);
		for (size_t j = 0; j < col_num_; j++)
			row[j] = data_[j][i];
		return row;
	}

	size_t RowSize() const
	{
		return row_num_;
	}

	size_t ColumnSize() const
	{
		return col_num_;
	}

	typename std::vector<std::vector<double>>::iterator begin()
	{
		return data_.begin();
	}

	typename std::vector<std::vector<double>>::const_iterator begin() const
	{
		return data_.begin();
	}

	typename std::vector<std::vector<double>>::iterator end()
	{
		return data_.end();
	}

	typename std::vector<std::vector<double>>::const_iterator end() const
	{
		return data_.end();
	}

	std::vector<double>& operator[](const size_t& index)
	{
		return data_[index];
	}

	const std::vector<double>& operator[](const size_t& index) const
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


Matrix getSubmatrix(const Matrix& M, const index_set& si)
{
	Matrix Res(M.RowSize(), 0);
	for (size_t i : si)
		Res.PushBackColumn(M[i]);
	return Res;
}


vector getSubvector(const vector& v, const index_set& si)
{
	vector res;
	for (size_t i : si)
		res.push_back(v[i]);
	return res;
}


index_set getNotExtendedBaseJ(const index_set& jbe, const size_t n)
{
	index_set res;
	size_t cur = 0;
	for (const size_t j : jbe)
	{
		while (cur < j)
			res.insert(cur++);
		cur++;
	}
	while (cur < n)
		res.insert(cur++);
	return res;
}


Matrix getMatrixH(const Matrix& A, const Matrix& D, const index_set& jbe)
{
	const size_t size1 = jbe.size(), size2 = A.RowSize();
	Matrix H(size1 + size2, size1 + size2);
	size_t i = 0, j = 0;
	for (size_t jbi : jbe)
	{
		j = 0;
		for (size_t jbj : jbe)
		{
			H[j][i] = D[jbj][jbi];
			j++;
		}
		i++;
	}
	for (i = 0; i < size2; i++)
	{
		j = 0;
		for (size_t jbj : jbe)
		{
			H[size1 + i][j] = A[jbj][i];
			H[j][size1 + i] = A[jbj][i];
			j++;
		}
	}
	for (i = 0; i < size2; i++)
		for (j = 0; j < size2; j++)
			H[size1 + j][size1 + i] = 0;
	return H;
}


vector getVectorL(const Matrix& A, const Matrix& D, const index_set& jbe, size_t j0)
{
	vector l(A.ColumnSize(), 0);
	l[j0] = 1;
	Matrix H = getMatrixH(A, D, jbe);
	vector bb;
	for (size_t j : jbe)
		bb.push_back(D[j0][j]);
	for (const double& el : A[j0])
		bb.push_back(el);
	vector xx = -(H.Inverse() * bb);
	size_t i = 0;
	for (size_t j : jbe)
		l[j] = xx[i++];
	return l;
}


bool Solve(Matrix& A, vector& b, vector& c, Matrix& D, vector& x, index_set& jb, index_set& jbe)
{
	size_t n = A.ColumnSize();
	index_set jben = getNotExtendedBaseJ(jbe, n);
	Matrix B(0, 0);
	bool is_eval_B = false;

	while (true)
	{
		if (!is_eval_B)
		{
			B = getSubmatrix(A, jb).Inverse();
			is_eval_B = true;
		}
		vector cx = c + D * x;
		vector ux = -getSubvector(cx, jb) * B;
		vector diff = ux * A + cx;

		size_t j0 = n;
		double min_diff = INF;
		for (const size_t j : jben)
			if (diff[j] < min_diff)
			{
				min_diff = diff[j];
				j0 = j;
			}
		if (min_diff >= -EPS)
			return true;

		while (true)
		{
			vector l = getVectorL(A, D, jbe, j0);
			size_t jp = j0;
			double delta = l * D * l, min_teta = INF;
			if (delta > EPS)
				min_teta = fabs(diff[j0]) / delta;
			for (size_t j : jbe)
			{
				if (l[j] < -EPS && -x[j] / l[j] < min_teta)
				{
					min_teta = -x[j] / l[j];
					jp = j;
				}
			}
			if (min_teta == INF)
				return false;

			x += min_teta * l;

			if (jp == j0)
			{
				jbe.insert(j0);
				jben.erase(j0);
				break;
			}
			else if (jb.find(jp) == jb.end() && jbe.find(jp) != jbe.end())
			{
				jbe.erase(jp);
				jben.insert(jp);
				diff[j0] += min_teta * delta;
			}
			else
			{
				size_t s = std::distance(jb.begin(), jb.find(jp));
				bool is_complete = false;
				if (jb.size() != jbe.size())
				{
					std::vector<size_t> diff_v;
					std::set_difference(jbe.begin(), jbe.end(), jb.begin(), jb.end(), std::back_inserter(diff_v));
					if (!is_eval_B)
					{
						B = getSubmatrix(A, jb).Inverse();
						is_eval_B = true;
					}
					vector s_row;
					for (size_t i = 0; i < B.ColumnSize(); i++)
						s_row.push_back(B[i][s]);
					for (size_t j : diff_v)
						if (fabs(s_row * A[j]) > EPS)
						{
							jb.erase(jp);
							jb.insert(j);
							jbe.erase(jp);
							jben.insert(jp);
							diff[j0] += min_teta * delta;
							is_eval_B = false;
							is_complete = true;
							break;
						}
				}
				if (!is_complete)
				{
					jb.erase(jp);
					jb.insert(j0);
					jbe.erase(jp);
					jbe.insert(j0);
					jben.erase(j0);
					jben.insert(jp);
					is_eval_B = false;
					break;
				}
			}
		}
	}

	return true;
}


int main()
{
	std::ios_base::sync_with_stdio(false);
	std::cin.tie(0);

	size_t m, n, k;
	std::cin >> m >> n >> k;

	Matrix A(m, n);
	std::cin >> A;

	vector b(m);
	std::cin >> b;

	vector c(n);
	std::cin >> c;

	Matrix D(n, n);
	std::cin >> D;

	vector x(n);
	std::cin >> x;

	index_set jb;
	for (size_t i = 0, j; i < m; i++)
	{
		std::cin >> j;
		jb.insert(j - 1);
	}

	index_set jbe;
	for (size_t i = 0, j; i < m; i++)
	{
		std::cin >> j;
		jbe.insert(j - 1);
	}

	try
	{
		bool is_bounded = Solve(A, b, c, D, x, jb, jbe);

		if (is_bounded)
		{
			std::cout << "Bounded\n";
			std::cout << std::fixed << std::setprecision(9);
			std::cout << x;
		}
		else
		{
			std::cout << "Unbounded";
		}
	}
	catch (std::exception & e)
	{
		std::cout << e.what();
	}

	return 0;
}
