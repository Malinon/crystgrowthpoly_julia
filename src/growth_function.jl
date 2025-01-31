
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

end # module gf