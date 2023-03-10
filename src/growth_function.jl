
module gf
using Printf
struct growth_function
    polynomials_coefficients::Vector{Integer}
    cell_dimmension::Integer
end

function get_description(self::growth_function)
    fun = @sprintf("f_%d", self.cell_dimmension)
    print(fun, "=", self.polynomials_coefficients)
end
end
# print(fun, "(n)=", polynomials_coefficients[0][0] reduce( (x,y) -> x*" + "*y, polynomials_coefficients[0][i] * "x^" * str(i) for i in 1:length(polynomials_coefficients[0])))