from utils import expect

from collections import OrderedDict
import pathlib, re

#
# Global hardcoded data
#

# Templates: maps piece name to generic file text
FILE_TEMPLATES = {
    "cxx_func_impl": lambda phys, sub, gen_code:
"""
#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/{physics}/{physics}_functions.hpp"
#include "physics/{physics}/{phyiscs}_functions_f90.hpp"

#include "{physics}_unit_tests_common.hpp"

namespace scream {
namespace {physics} {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::{test_data_struct} {

{gen_code}

}; // struct {test_data_struct}

} // namespace unit_test
} // namespace {physics}
} // namespace scream

namespace {

TEST_CASE({sub}_bfb, "[{physics}_functions]")
{
  using TestStruct = scream::{physics}::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::{test_data_struct};

  TestStruct::run_bfb();
}

} // empty namespace
""".format(physics=phys, test_data_struct=get_data_test_struct_name(sub), gen_code=gen_code),

    "cxx_bfb_unit_impl": lambda phys, sub, gen_code:
"""
#ifndef {phys_upper}_{sub_upper}_IMPL_HPP
#define {phys_upper}_{sub_upper}_IMPL_HPP

#include "{physics}_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace {physics} {

/*
 * Implementation of {physics} {sub}. Clients should NOT
 * #include this file, but include {physics}_functions.hpp instead.
 */

{gen_code}

} // namespace p3
} // namespace scream

#endif
""".format(physics=phys, test_data_struct=get_data_test_struct_name(sub), gen_code=gen_code, phys_upper=phys.upper(), sub_upper=sub.upper()),

    "cxx_eti": lambda phys, sub, gen_code: gen_code
}

# piece map. maps the name of a piece of boilerplate that needs to be genereate to:
#   (filepath (relative to cxx_root))
FILEPATH, FILECREATE, INSERT_REGEX, ID_SELF_BEGIN_REGEX, ID_SELF_END_REGEX = range(5)
PIECES = {
    "f90_c2f_bind"  : (
        lambda phys, sub, gb: "{}_iso_c.f90".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "f90_c2f_bind"),
        lambda phys, sub, gb: re.compile(r"^\s*end\s+module\s{}_iso_c".format(phys)), # put at end of module
        lambda phys, sub, gb: get_subroutine_begin_regex(sub + "_c"),
        lambda phys, sub, gb: get_subroutine_end_regex(sub + "_c"),
    ),

    "f90_f2c_bind"  : (
        lambda phys, sub, gb: "{}_iso_f.f90".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "f90_f2c_bind"),
        lambda phys, sub, gb: re.compile(r"^\s*end\s+module\s{}_iso_f".format(phys)), # put at end of module
        lambda phys, sub, gb: get_subroutine_begin_regex(sub + "_f"),
        lambda phys, sub, gb: get_subroutine_end_regex(sub + "_f"),
    ),

    "cxx_c2f_bind_decl"  : (
        lambda phys, sub, gb: "{}_functions_f90.cpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_bind_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(comment='extern "C" : end _c decls'), # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_c"),
        lambda phys, sub, gb: re.compile(r".*;\s*$"),
    ),

    "cxx_c2f_glue_decl"  : (
        lambda phys, sub, gb: "{}_functions_f90.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_glue_decl"),
        lambda phys, sub, gb: re.compile(r'^\s*extern\s+"C"'), # put before _f decls
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub),
        lambda phys, sub, gb: re.compile(r".*"),
    ),

    "cxx_c2f_glue_impl"  : (
        lambda phys, sub, gb: "{}_functions_f90.cpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_glue_impl"),
        lambda phys, sub, gb: re.compile(r"^\s*// end _c impls"), # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub),
        lambda phys, sub, gb: re.compile(r"^[}]\s*$"),
    ),

    "cxx_c2f_data"  : (
        lambda phys, sub, gb: "{}_functions_f90.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_c2f_data"),
        lambda phys, sub, gb: re.compile(r"^\s*// Glue functions to call fortran"), # reqs special comment
        lambda phys, sub, gb: get_cxx_struct_begin_regex(get_data_struct_name(sub)),
        lambda phys, sub, gb: re.compile(r".*(namespace|struct|// Glue)"),
    ),

    "cxx_f2c_bind_decl"  : (
        lambda phys, sub, gb: "{}_functions_f90.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_f2c_bind_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(comment="end _f function decls"), # reqs special comment
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_f"),
        lambda phys, sub, gb: re.compile(r".*;\s*$"),
    ),

    "cxx_f2c_bind_impl"  : (
        lambda phys, sub, gb: "{}_functions_f90.cpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_f2c_bind_impl"),
        lambda phys, sub, gb: get_namespace_close_regex(phys),
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub + "_f"),
        lambda phys, sub, gb: get_cxx_close_block_regex(comment="end {}".format(sub + "_f")), # reqs special comment
    ),

    "cxx_func_decl"  : (
        lambda phys, sub, gb: "{}_functions.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_func_hpp"),
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True, comment="struct Functions"),
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub),
        lambda phys, sub, gb: re.compile(r".*;\s*$"),
    ),

    "cxx_func_impl" : (
        lambda phys, sub, gb: "{}_{}_impl.hpp".format(phys, sub),
        lambda phys, sub, gb: create_template(phys, sub, gb, "cxx_func_impl"),
        lambda phys, sub, gb: get_namespace_close_regex(phys),
        lambda phys, sub, gb: get_cxx_function_begin_regex(sub),
        lambda phys, sub, gb: re.compile(r".*{\s*$"),
    ),

    "cxx_bfb_unit_decl" : (
        lambda phys, sub, gb: "tests/{}_unit_tests_common.hpp".format(phys),
        lambda phys, sub, gb: expect_exists(phys, sub, gb, "cxx_bfb_unit_decl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True),
        lambda phys, sub, gb: get_cxx_struct_begin_regex(get_data_test_struct_name(sub)),
        lambda phys, sub, gb: re.compile(r".*;\s*$"),
    ),

    "cxx_bfb_unit_impl"  : (
        lambda phys, sub, gb: "tests/{}_{}_tests.cpp".format(phys, sub),
        lambda phys, sub, gb: create_template(phys, sub, gb, "cxx_bfb_unit_impl"),
        lambda phys, sub, gb: get_cxx_close_block_regex(semicolon=True, comment="struct {}".format(get_data_test_struct_name(sub))),
        lambda phys, sub, gb: get_cxx_function_begin_regex("run_bfb"),
        lambda phys, sub, gb: get_cxx_close_block_regex(comment="run_bfb"),
    ),

    "cxx_eti" : (
        lambda phys, sub, gb: "{}_{}.cpp".format(phys, sub),
        lambda phys, sub, gb: create_template(phys, sub, gb, "cxx_eti"),
        lambda phys, sub, gb: re.compile(".*"),
        lambda phys, sub, gb: re.compile(".*"),
        lambda phys, sub, gb: get_namespace_close_regex("scream"),
    ),
}

# physics map. maps the name of a physics packages containing the original fortran subroutines to:
#   (path-to-origin, path-to-cxx-src)
ORIGIN_FILE, CXX_ROOT, INIT_CODE = range(3)
PHYSICS = {
    "p3"   : (
        "components/cam/src/physics/cam/micro_p3.F90",
        "components/scream/src/physics/p3",
        "p3_init();"
    ),
    "shoc" : (
        "components/cam/src/physics/cam/shoc.F90",
        "components/scream/src/physics/shoc",
        "shoc_init(d.nlev, true);"
    ),
}

#
# Free functions
#

###############################################################################
def get_subroutine_begin_regex(name):
###############################################################################
    subroutine_begin_regex_str = r"^\s*subroutine\s+{}\s*[(]".format(name)
    return re.compile(subroutine_begin_regex_str)

###############################################################################
def get_subroutine_end_regex(name):
###############################################################################
    subroutine_end_regex_str = r"^\s*end\s+subroutine\s+{}\s*$".format(name)
    return re.compile(subroutine_end_regex_str)

###############################################################################
def get_cxx_function_begin_regex(name):
###############################################################################
    function_begin_regex_str = r"^\s*void\s+{}\s*[(]".format(name)
    return re.compile(function_begin_regex_str)

###############################################################################
def get_cxx_close_block_regex(semicolon=False, comment=None):
###############################################################################
    semicolon_regex_str = r"\s*;" if semicolon else ""
    comment_regex_str   = r"\s*//\s*{}".format(comment) if comment else ""
    close_block_regex_str = re.compile(r"^\s*}{}{}\s*$".format(semicolon_regex_str, comment_regex_str))
    return re.compile(close_block_regex_str)

###############################################################################
def get_namespace_close_regex(namespace):
###############################################################################
    return get_cxx_close_block_regex(comment=r"namespace\s+{}".format(namespace))

###############################################################################
def get_cxx_struct_begin_regex(struct):
###############################################################################
    struct_regex_str = r"^\s*struct\s+{}([\W]|$)".format(struct)
    return re.compile(struct_regex_str)

###############################################################################
def get_data_struct_name(sub):
###############################################################################
    return "".join([item.capitalize() for item in sub.split("_")])

###############################################################################
def get_data_test_struct_name(sub):
###############################################################################
    return "Test{}".format(get_data_struct_name(sub))

###############################################################################
def get_supported_pieces():
###############################################################################
    return PIECES.keys()

###############################################################################
def get_supported_physics():
###############################################################################
    return PHYSICS.keys()

###############################################################################
def get_piece_data(physics, sub, piece_name, piece_data):
###############################################################################
    return PIECES[piece_name][piece_data](physics, sub)

###############################################################################
def get_physics_data(physics_name, physics_data):
###############################################################################
    return PHYSICS[physics_name][physics_data]

###############################################################################
def expect_exists(physics, sub, piece, gb):
###############################################################################
    filepath = gb.get_path_for_piece_file(physics, sub, piece, gb)
    expect(filepath.exists(), "For generating {}'s {} for phyiscs {}, expected file {} to already exist".\
           format(sub, piece, physics, filepath))
    return False # File was not created

###############################################################################
def create_template(physics, sub, piece, gb):
###############################################################################
    """
    Create a file based on a template if it doesn't exist. Return True if a file was created.
    """
    filepath = gb.get_path_for_piece_file(physics, sub, piece)
    if not filepath.exists():
        expect(piece in FILE_TEMPLATES,
               "{} does not exist and there is no template for generating files for piece {}".format(filepath, piece))

        gen_code = getattr(gb, "gen_{}".format(piece))(physics, sub)
        contents = FILE_TEMPLATES[piece](physics, sub, gen_code)
        if gb.dry_run():
            print("Would create file {} with contents:\n\n{}".format(filepath, contents))
        else:
            with filepath.open("w") as fd:
                fd.write(contents)

        return True
    else:
        return False

###############################################################################
def remove_comments_and_ws(contents):
###############################################################################
    """
    >>> teststr = '''
    ... module mymod
    ...   subroutine foo(a, b, &{0}
    ...                c, d, e,&
    ... !bad{0}
    ... &f)
    ...
    ...     real, intent(in) :: a, b, & !go
    ...                 c, d, e, f{0}
    ...
    ...   ! hi
    ...   !hi ! there{0}
    ... !hi ! there
    ...   end subroutine foo{0}
    ... end module mymod
    ... '''.format(" ")
    >>> print(remove_comments_and_ws(teststr))
    module mymod
    subroutine foo(a, b, &
    c, d, e,&
    &f)
    real, intent(in) :: a, b, &
    c, d, e, f
    end subroutine foo
    end module mymod
    """
    new_lines = []
    comment_regex = re.compile(r"^([^!]*)")
    for line in contents.splitlines():
        m = comment_regex.match(line)
        if m is not None:
            line = m.groups()[0].strip()
        else:
            line = line.strip()

        if line != "":
            new_lines.append(line)

    return "\n".join(new_lines)

###############################################################################
def resolve_line_continuations(contents):
###############################################################################
    """
    >>> teststr = '''
    ... module mymod
    ...   subroutine foo(a, b, &
    ...                c, d, e,&{0}
    ... !bad
    ... &f){0}
    ...
    ...     real, intent(in) :: a, b, & !go{0}
    ...                 c, d, e, f
    ...
    ...   ! hi
    ...   !hi ! there
    ... !hi ! there{0}
    ...   end subroutine foo
    ... end module mymod{0}
    ... '''.format("  ")
    >>> print(resolve_line_continuations(teststr))
    module mymod
    subroutine foo(a, b, c, d, e,f)
    real, intent(in) :: a, b, c, d, e, f
    end subroutine foo
    end module mymod
    """
    # Must remove comments and whitespace for the alg below to work
    contents = remove_comments_and_ws(contents)

    new_lines = []
    contination = False
    for line in contents.splitlines():
        line = line.lstrip("&") # always remove leading &, they are optional
        if line != "":
            if contination:
                new_lines[-1] += line.rstrip("&")
            else:
                new_lines.append(line.rstrip("&"))

            contination = line.endswith("&")

    return "\n".join(new_lines)

###############################################################################
def normalize_f90(contents):
###############################################################################
    # We do not currently attempt to preprocess the contents
    return resolve_line_continuations(contents).lower()

ARG_NAME, ARG_TYPE, ARG_INTENT, ARG_DIMS = range(4)
###############################################################################
def parse_f90_args(line):
###############################################################################
    """
    Given a line of fortran code declaring an argument[s], return [(argname, argtype, intent, dims)]

    >>> parse_f90_args('integer, intent(in) :: kts, kte, kbot')
    [('kts', 'integer', 'in', None), ('kte', 'integer', 'in', None), ('kbot', 'integer', 'in', None)]
    >>> parse_f90_args('real(rtype),intent(inout ), dimension(kts:kte) :: pres,dpres,  dz ')
    [('pres', 'real', 'inout', ('kts:kte',)), ('dpres', 'real', 'inout', ('kts:kte',)), ('dz', 'real', 'inout', ('kts:kte',))]
    >>> parse_f90_args('logical (btype), intent( in) ::do_predict_nc')
    [('do_predict_nc', 'logical', 'in', None)]
    >>> parse_f90_args('real(rtype),intent(inout), dimension( kts:kte, its: ite) :: dz')
    [('dz', 'real', 'inout', ('kts:kte', 'its:ite'))]
    >>> parse_f90_args('real(rtype),intent(inout), dimension(3) :: dz')
    [('dz', 'real', 'inout', ('3',))]
    >>> parse_f90_args('real(rtype),intent(inout), dimension(3,4) :: dz')
    [('dz', 'real', 'inout', ('3', '4'))]
    >>> parse_f90_args('real(rtype), dimension(3,4),intent(inout) :: dz')
    [('dz', 'real', 'inout', ('3', '4'))]
    """
    expect(line.count("::") == 1, "Expected line format 'type-info :: names' for: {}".format(line))
    metadata, names_str = line.split("::")
    names = [name.strip() for name in names_str.split(",")]
    metadata_raw = [item.strip() for item in metadata.strip().split(",")]
    metadata = []
    balanced = True
    for metadatum_raw in metadata_raw:
        if balanced:
            metadata.append(metadatum_raw)
        else:
            metadata[-1] += ",{}".format(metadatum_raw)

        balanced = metadata[-1].count("(") == metadata[-1].count(")")

    argtype = metadata[0].split("(")[0].strip()
    intent, dims = None, None
    for metadatum in metadata:
        if metadatum.startswith("intent"):
            expect(intent is None, "Multiple intents in line: {}".format(line))
            intent = metadatum.split("(")[-1].rstrip(")").strip()
        elif metadatum.startswith("dimension"):
            expect(dims is None, "Multiple dimensions in line: {}".format(line))
            dims_raw = metadatum.split("(")[-1].rstrip(")").strip()
            dims = tuple(item.replace(" ", "") for item in dims_raw.split(","))

    return [(name, argtype, intent, dims) for name in names]

###############################################################################
def parse_origin(contents, subs):
###############################################################################
    """
    Returns a map of subname->[(argname, argtype, intent, dims)]

    >>> teststr = '''
    ...
    ...   SUBROUTINE p3_get_tables(mu_r_user, revap_user, vn_user, vm_user)
    ...     ! This can be called after p3_init_b.
    ...     implicit none
    ...     real(rtype), dimension(150), intent(out) :: mu_r_user
    ...     real(rtype), dimension(300,10), intent(out) :: vn_user, vm_user, revap_user
    ...     mu_r_user(:) = mu_r_table(:)
    ...     revap_user(:,:) = revap_table(:,:)
    ...     vn_user(:,:) = vn_table(:,:)
    ...     vm_user(:,:) = vm_table(:,:)
    ...
    ...    return
    ...
    ...   end SUBROUTINE p3_get_tables
    ...
    ...   subroutine p3_set_tables(mu_r_user, revap_user, vn_user, vm_user)
    ...     ! This can be called instead of p3_init_b.
    ...     implicit none
    ...     real(rtype), dimension(150), intent(in) :: mu_r_user
    ...     real(rtype), dimension(300,10), intent(in) :: vn_user, vm_user, revap_user
    ...     mu_r_table(:) = mu_r_user(:)
    ...     revap_table(:,:) = revap_user(:,:)
    ...     vn_table(:,:) = vn_user(:,:)
    ...     vm_table(:,:) = vm_user(:,:)
    ...
    ...    return
    ...
    ...   END SUBROUTINE p3_set_tables
    ...
    ...   SUBROUTINE p3_init_b()
    ...     implicit none
    ...     integer                      :: i,ii,jj,kk
    ...     real(rtype)                         :: lamr,mu_r,dm,dum1,dum2,dum3,dum4,dum5,  &
    ...          dd,amg,vt,dia
    ...
    ...     ! AaronDonahue: Switching to table ver 4 means switching to a constand mu_r,
    ...     ! so this section is commented out.
    ...     do i = 1,150
    ...   END SUBROUTINE p3_init_b
    ... '''
    >>> sorted(parse_origin(teststr, ["p3_get_tables", "p3_init_b"]).items())
    [('p3_get_tables', [('mu_r_user', 'real', 'out', ('150',)), ('vn_user', 'real', 'out', ('300', '10')), ('vm_user', 'real', 'out', ('300', '10')), ('revap_user', 'real', 'out', ('300', '10'))]), ('p3_init_b', [])]
    """
    begin_regexes = [get_subroutine_begin_regex(sub) for sub in subs]
    arg_decl_regex = re.compile(r"^.+intent\s*[(]\s*(in|out|inout)\s*[)]")

    contents = normalize_f90(contents)

    db = {}
    active_sub = None
    arg_decls = []
    for line in contents.splitlines():
        begin_match = None
        for sub, begin_regex in zip(subs, begin_regexes):
            begin_match = begin_regex.match(line)
            if begin_match is not None:
                expect(active_sub is None, "subroutine {} was still active when {} began".format(active_sub, sub))
                active_sub = sub

        if active_sub:
            decl_match = arg_decl_regex.match(line)
            if decl_match is not None:
                arg_decls.extend(parse_f90_args(line))

            end_regex = get_subroutine_end_regex(active_sub)
            end_match = end_regex.match(line)
            if end_match is not None:
                expect(active_sub not in db, "Found multiple matches for {}".format(active_sub))
                db[active_sub] = arg_decls
                active_sub = None
                arg_decls = []

    return db

C_TYPE_MAP = {"real" : "c_real", "integer" : "c_int", "logical" : "c_bool"}
###############################################################################
def gen_arg_f90_decl(argtype, intent, dims, names):
###############################################################################
    """
    Generate f90 argument declaration based on the input data

    >>> gen_arg_f90_decl("real", "in", ("10", "150"), ["foo", "bar"])
    'real(kind=c_real) , intent(in), dimension(10, 150) :: foo, bar'
    >>> gen_arg_f90_decl("real", "out", ("10", "150"), ["foo", "bar"])
    'real(kind=c_real) , intent(out), dimension(10, 150) :: foo, bar'
    >>> gen_arg_f90_decl("logical", "in", None, ["biz", "baz"])
    'logical(kind=c_bool) , value, intent(in) :: biz, baz'
    >>> gen_arg_f90_decl("integer", "inout", None, ["barg"])
    'integer(kind=c_int) , intent(inout) :: barg'
    >>> gen_arg_f90_decl("integer", "out", None, ["barg"])
    'integer(kind=c_int) , intent(out) :: barg'
    """
    expect(argtype in C_TYPE_MAP, "Unrecognized argtype for C_TYPE_MAP: {}".format(argtype))
    c_type = C_TYPE_MAP[argtype]
    value  = ", value" if dims is None and intent == "in" else ""
    intent_s = ", intent({})".format(intent)
    dimension_s = ", dimension({})".format(", ".join(dims)) if dims is not None else ""
    names_s = ", ".join(names)
    return "{argtype}(kind={c_type}) {value}{intent}{dimension} :: {names}".\
        format(argtype=argtype, c_type=c_type, value=value, intent=intent_s, dimension=dimension_s, names=names_s)

CXX_TYPE_MAP = {"real" : "Real", "integer" : "Int", "logical" : "bool"}
###############################################################################
def get_cxx_type(arg_datum):
###############################################################################
    """
    Based on arg datum, give c++ bridge type

    >>> get_cxx_type(("foo", "real", "in", ("100",)))
    'Real*'
    >>> get_cxx_type(("foo", "real", "in", None))
    'Real'
    >>> get_cxx_type(("foo", "real", "inout", ("100",)))
    'Real*'
    >>> get_cxx_type(("foo", "real", "inout", None))
    'Real*'
    >>> get_cxx_type(("foo", "real", "out", ("100",)))
    'Real*'
    >>> get_cxx_type(("foo", "real", "out", None))
    'Real*'
    >>> get_cxx_type(("foo", "integer", "inout", None))
    'Int*'
    """
    is_ptr = arg_datum[ARG_DIMS] is not None or arg_datum[ARG_INTENT] != "in"
    arg_type = arg_datum[ARG_TYPE]
    expect(arg_type in CXX_TYPE_MAP, "Unrecognized argtype for CXX_TYPE_MAP: {}".format(arg_type))
    arg_cxx_type = CXX_TYPE_MAP[arg_type]
    return "{}{}".format(arg_cxx_type, "*" if is_ptr else "")

KOKKOS_TYPE_MAP = {"real" : "Spack", "integer" : "Int", "logical" : "bool"}
###############################################################################
def get_kokkos_type(arg_datum):
###############################################################################
    """
    Based on arg datum, give c++ kokkos type

    Note: We can only guess at the correct types, especially whether an argument
    should be packed data or not!

    >>> get_kokkos_type(("foo", "real", "in", ("100",)))
    'const uview_1d<const Spack>&'
    >>> get_kokkos_type(("foo", "real", "in", None))
    'const Spack&'
    >>> get_kokkos_type(("foo", "real", "inout", ("100",)))
    'const uview_1d<Spack>&'
    >>> get_kokkos_type(("foo", "real", "inout", None))
    'Spack&'
    >>> get_kokkos_type(("foo", "real", "out", ("100",)))
    'const uview_1d<Spack>&'
    >>> get_kokkos_type(("foo", "real", "out", None))
    'Spack&'
    >>> get_kokkos_type(("foo", "integer", "inout", None))
    'Int&'
    """
    is_const  = arg_datum[ARG_INTENT] == "in"
    is_view   = arg_datum[ARG_DIMS] is not None
    base_type = "{}{}".format("const " if is_const else "", KOKKOS_TYPE_MAP[arg_datum[ARG_TYPE]])

    # We assume 1d even if the f90 array is 2d since we assume c++ will spawn a kernel
    # over one of the dimensions
    return "const uview_1d<{}>&".format(base_type) if is_view else "{}&".format(base_type)

###############################################################################
def gen_arg_cxx_decls(arg_data, kokkos=False):
###############################################################################
    """
    Get all arg decls for a set of arg data

    >>> gen_arg_cxx_decls([("foo", "real", "in", ("100",)), ("bar", "real", "in", None)])
    ['Real* foo', 'Real bar']
    >>> gen_arg_cxx_decls([("foo", "real", "in", ("100",)), ("bar", "real", "in", None)], kokkos=True)
    ['const uview_1d<const Spack>& foo', 'const Spack& bar']
    """
    arg_names    = [item[ARG_NAME] for item in arg_data]
    get_type     = get_kokkos_type if kokkos else get_cxx_type
    arg_types    = [get_type(item) for item in arg_data]
    arg_sig_list = ["{} {}".format(arg_type, arg_name) for arg_name, arg_type in zip(arg_names, arg_types)]
    return arg_sig_list

###############################################################################
def needs_transpose(arg_data):
###############################################################################
    """
    Based on data, does this sub need to have its data transposed when going between languages?

    >>> needs_transpose([("foo", "real", "in", ("100",)), ("bar", "real", "in", None)])
    False
    >>> needs_transpose([("foo", "real", "in", ("100","50")), ("bar", "real", "in", None)])
    True
    >>> needs_transpose([("foo", "real", "in", None), ("bar", "real", "in", None)])
    False
    """
    for arg_datum in arg_data:
        arg_dims = arg_datum[ARG_DIMS]
        if arg_dims is not None and len(arg_dims) > 1:
            return True

    return False

###############################################################################
def gen_cxx_data_args(arg_data):
###############################################################################
    """
    Based on data, generate unpacking of Data struct args

    >>> gen_cxx_data_args([("foo", "real", "in", ("100",)), ("bar", "real", "in", None), ("baz", "integer", "inout", None), ("bag", "integer", "inout", ("100",))])
    ['d.foo', 'd.bar', '&d.baz', 'd.bag']
    """
    args_needs_ptr = [item[ARG_DIMS] is None and item[ARG_INTENT] != "in" for item in arg_data]
    arg_names      = [item[ARG_NAME] for item in arg_data]
    args = ["{}d.{}".format("&" if need_ptr else "", arg_name) for arg_name, need_ptr in zip(arg_names, args_needs_ptr)]
    return args

###############################################################################
def gen_arg_f90_decls(arg_data):
###############################################################################
    r"""
    Generate f90 argument declarations, will attempt to group these together if possible.

    >>> arg_data = [
    ... ("foo1", "real", "in", ("100",)),
    ... ("foo2", "real", "in", ("100",)),
    ... ("bar1", "real", "in", ("100","50")),
    ... ("bar2", "real", "in", ("100","50")),
    ... ("baz", "real", "inout", ("100",)),
    ... ("bag", "integer", "in", ("100",)),
    ... ("bab1", "integer", "out", None),
    ... ("bab2", "integer", "out", None),
    ... ("val", "logical", "in", None),
    ... ]
    >>> print("\n".join(gen_arg_f90_decls(arg_data)))
    real(kind=c_real) , intent(in), dimension(100) :: foo1, foo2
    real(kind=c_real) , intent(in), dimension(100, 50) :: bar1, bar2
    real(kind=c_real) , intent(inout), dimension(100) :: baz
    integer(kind=c_int) , intent(in), dimension(100) :: bag
    integer(kind=c_int) , intent(out) :: bab1, bab2
    logical(kind=c_bool) , value, intent(in) :: val
    """
    metadata = OrderedDict()
    for name, argtype, intent, dims in arg_data:
        metatuple = (argtype, intent, dims)
        metadata.setdefault(metatuple, []).append(name)

    result = []
    for metatuple, names in metadata.items():
        result.append(gen_arg_f90_decl(*metatuple, names))

    return result

###############################################################################
def has_arrays(arg_data):
###############################################################################
    for _, _, _, dims in arg_data:
        if dims is not None:
            return True

    return False

###############################################################################
def gen_struct_members(arg_data):
###############################################################################
    r"""
    Gen cxx code for data struct members

    >>> arg_data = [
    ... ("foo1", "real", "in", ("100",)),
    ... ("foo2", "real", "in", ("100",)),
    ... ("bar1", "real", "in", ("100","50")),
    ... ("bar2", "real", "in", ("100","50")),
    ... ("gag", "real", "in", None),
    ... ("baz", "real", "inout", ("100",)),
    ... ("bag", "integer", "in", ("100",)),
    ... ("bab1", "integer", "out", None),
    ... ("bab2", "integer", "out", None),
    ... ("val", "logical", "in", None),
    ... ]
    >>> print("\n".join(gen_struct_members(arg_data)))
    // Inputs
    Real *foo1, *foo2, *bar1, *bar2;
    Real gag;
    Int *bag;
    bool val;
    <BLANKLINE>
    // Inputs/Outputs
    Real *baz;
    <BLANKLINE>
    // Outputs
    Int bab1, bab2;
    <BLANKLINE>
    """
    metadata = {} # intent -> (type, is_ptr) -> names
    for name, argtype, intent, dims in arg_data:
        metadata.setdefault(intent, OrderedDict()).setdefault((argtype, dims is not None), []).append(name)

    intent_order = ( ("in", "Inputs"), ("inout", "Inputs/Outputs"), ("out", "Outputs") )
    result = []
    for intent, comment in intent_order:
        if intent in metadata:
            result.append("// {}".format(comment))
            type_map = metadata[intent]
            for type_info, names in type_map.items():
                type_name, is_ptr = type_info
                decl_str = CXX_TYPE_MAP[type_name]
                decl_str += " {};".format(", ".join(["{}{}".format("*" if is_ptr else "", name) for name in names]))
                result.append(decl_str)

            result.append("")

    return result

###############################################################################
def group_data(arg_data, filter_out_intent=None):
###############################################################################
    """
    Given data, return (i_dim, k_dim, j_dim, [scalars], [ik_reals], [ij_reals], [i_reals], [i_ints])

    >>> arg_data = [
    ... ("foo1", "real", "in", ("shcol",)),
    ... ("foo2", "real", "in", ("shcol",)),
    ... ("bar1", "real", "in", ("shcol","nlev")),
    ... ("bar2", "real", "in", ("shcol","nlev")),
    ... ("bak1", "real", "in", ("shcol","nlevi")),
    ... ("bak2", "real", "in", ("shcol","nlevi")),
    ... ("gag", "real", "in", None),
    ... ("baz", "real", "inout", ("shcol",)),
    ... ("bag", "integer", "in", ("shcol",)),
    ... ("bab1", "integer", "out", None),
    ... ("bab2", "integer", "out", None),
    ... ("val", "logical", "in", None),
    ... ("shcol", "integer", "in", None),
    ... ("nlev", "integer", "in", None),
    ... ("nlevi", "integer", "in", None),
    ... ]
    >>> group_data(arg_data)
    ('shcol', 'nlev', 'nlevi', [('gag', 'Real'), ('bab1', 'Int'), ('bab2', 'Int'), ('val', 'bool')], ['bar1', 'bar2'], ['bak1', 'bak2'], ['foo1', 'foo2', 'baz'], ['bag'])
    """
    i_ints   = []
    i_reals  = []
    ik_reals = []
    ij_reals = [] # j = alternate k dim
    scalars  = []

    i_dim = None
    k_dim = None
    j_dim = None

    for name, argtype, _, dims in arg_data:
        if dims is not None:
            expect(len(dims) >= 1 and len(dims) <= 2,
                   "Only 1d and 2d data is support, {} has too many dims".format(name))

            if i_dim is None:
                i_dim = dims[0]
            else:
                expect(i_dim == dims[0],
                       "Ambiguous dimension for {}, expected {}, got{}".format(name, i_dim, dims[0]))

            if len(dims) == 2:
                if k_dim is None:
                    k_dim = dims[1]
                elif k_dim == dims[1]:
                    pass
                elif j_dim is None:
                    j_dim = dims[1]
                elif j_dim == dims[1]:
                    pass
                else:
                    expect(False, "Unable to identify dimension {} for arg {}".format(dims[1], name))

    for name, argtype, intent, dims in arg_data:
        if filter_out_intent is None or intent != filter_out_intent:
            if dims is None:
                if (name not in [i_dim, k_dim, j_dim]):
                    scalars.append( (name, CXX_TYPE_MAP[argtype]))
                else:
                    expect(argtype == "integer", "Expected dimension {} to be of type integer".format(name))

            elif argtype == "integer":
                expect(len(dims) == 1 and dims[0] == i_dim, "integer data {} has unsupported dims {}".format(name, dims))
                i_ints.append(name)

            elif len(dims) == 1:
                expect(dims[0] == i_dim and argtype in ["integer", "real"],
                       "1d real data {} has unsupported dims {}".format(name, dims))
                if argtype == "real":
                    i_reals.append(name)
                else:
                    i_ints.append(name)

            else:
                expect(len(dims) == 2 and dims[0] == i_dim and dims[1] in [k_dim, j_dim],
                       "2d real data {} has unsupported dims {}".format(name, dims))
                if dims[1] == k_dim:
                    ik_reals.append(name)
                else:
                    ij_reals.append(name)

    return i_dim, k_dim, j_dim, scalars, ik_reals, ij_reals, i_reals, i_ints

###############################################################################
def gen_struct_api(physics, struct_name, arg_data):
###############################################################################
    r"""
    Given data, generate code for data struct api

    >>> arg_data = [
    ... ("foo1", "real", "in", ("shcol",)),
    ... ("foo2", "real", "in", ("shcol",)),
    ... ("bar1", "real", "in", ("shcol","nlev")),
    ... ("bar2", "real", "in", ("shcol","nlev")),
    ... ("bak1", "real", "in", ("shcol","nlevi")),
    ... ("bak2", "real", "in", ("shcol","nlevi")),
    ... ("gag", "real", "in", None),
    ... ("baz", "real", "inout", ("shcol",)),
    ... ("bag", "integer", "in", ("shcol",)),
    ... ("bab1", "integer", "out", None),
    ... ("bab2", "integer", "out", None),
    ... ("val", "logical", "in", None),
    ... ("shcol", "integer", "in", None),
    ... ("nlev", "integer", "in", None),
    ... ("nlevi", "integer", "in", None),
    ... ]
    >>> print("\n".join(gen_struct_api("shoc", "DataSubName", arg_data)))
    DataSubName(Int shcol_, Int nlev_, Int nlevi_, Real gag_, Int bab1_, Int bab2_, bool val_) :
      PhysicsTestData(shcol_, nlev_, nlevi_, {&bar1, &bar2}, {&bak1, &bak2}, {&foo1, &foo2, &baz}, {&bag}), gag(gag_), bab1(bab1_), bab2(bab2_), val(val_) {}
    <BLANKLINE>
    SHOC_SCALARS(DataSubName, 3, 4, gag, bab1, bab2, val)
    """
    i_dim, k_dim, j_dim, scalars, ik_reals, ij_reals, i_reals, i_ints = group_data(arg_data)

    result = []
    dim_args = [(item, "Int") for item in [i_dim, k_dim, j_dim] if item is not None]
    cons_args = dim_args + scalars
    result.append("{struct_name}({cons_args}) :".\
                  format(struct_name=struct_name,
                         cons_args=", ".join(["{} {}_".format(argtype, name) for name, argtype in cons_args])))
    parent_call = "  PhysicsTestData({}".format(", ".join(["{}_".format(name) for name, _ in dim_args]))
    for item in (ik_reals, ij_reals, i_reals, i_ints):
        if len(item) > 0:
            parent_call += ", {{{}}}".format(", ".join(["&{}".format(name) for name in item]))

    parent_call += ")"
    if scalars:
        parent_call += ", {}".format(", ".join(["{0}({0}_)".format(name) for name, _ in scalars]))

    parent_call += " {}"
    result.append(parent_call)
    result.append("")

    if physics != "p3":
        if len(scalars) == 0:
            result.append("SHOC_NO_SCALAR({}, {})".format(struct_name, len(dim_args)))
        else:
            result.append("SHOC_SCALARS({}, {}, {}, {})".format(struct_name, len(dim_args), len(scalars),
                                                                ", ".join([name for name, _ in scalars])))
    else:
        expect(False, "p3 is not supported for now") # TODO

    return result

###############################################################################
def find_insertion(lines, insert_regex):
###############################################################################
    for idx, line in enumerate(lines):
        has_match = bool(insert_regex.match(line))
        if has_match:
            return idx

    return None

###############################################################################
def check_existing_piece(lines, begin_regex, end_regex):
###############################################################################
    begin_idx = None
    end_idx   = None

    for idx, line in enumerate(lines):
        begin_match = bool(begin_regex.match(line))
        end_match   = bool(end_regex.match(line))

        if begin_match:
            expect(begin_idx is None,
                   "Found multiple begin matches for pattern '{}' before end pattern '{}' was found".\
                   format(begin_regex.pattern, end_regex.pattern))

            if begin_idx is not None:
                begin_idx = idx

        if end_match and begin_idx is not None:
            end_idx = idx
            break

    if begin_idx is not None:
        expect(end_idx is not None,
               "Found no ending match for begin pattern '{}' starting on line {} and searching end pattern '{}'".\
               format(begin_regex.pattern, begin_idx, end_regex.pattern))

    return None if begin_idx is None else (begin_idx, end_idx)

#
# Main classes
#

###############################################################################
class GenBoiler(object):
###############################################################################

    ###########################################################################
    def __init__(self, subs, pieces, physics, overwrite, kernel, source_repo, target_repo, dry_run):
    ###########################################################################
        self._subs        = subs
        self._pieces      = pieces
        self._physics     = physics
        self._overwrite   = overwrite
        self._kernel      = kernel
        self._source_repo = pathlib.Path(source_repo).resolve()
        self._target_repo = pathlib.Path(target_repo).resolve()
        self._dry_run     = dry_run
        self._db          = {}

    ###########################################################################
    def _get_db(self, phys):
    ###########################################################################
        if phys in self._db:
            return self._db[phys]
        else:
            origin_file = self._source_repo / get_physics_data(phys, ORIGIN_FILE)
            expect(origin_file.exists(), "Missing origin file for physics {}: {}".format(phys, origin_file))
            db = parse_origin(origin_file.open().read(), self._subs)
            self._db[phys] = db
            return db

    ###########################################################################
    def _get_arg_data(self, phys, sub):
    ###########################################################################
        phys_db = self._get_db(phys)
        expect(sub in phys_db, "No data for subroutine {} in physics {}".format(sub, phys))
        return phys_db[sub]

    ###########################################################################
    def dry_run(self):
    ###########################################################################
        return self._dry_run

    ###############################################################################
    def get_path_for_piece_file(self, physics, sub, piece):
    ###############################################################################
        root_dir = pathlib.Path(get_physics_data(physics, CXX_ROOT))
        filepath = self._target_repo / root_dir / get_piece_data(physics, sub, piece, FILEPATH)
        return filepath

    ###########################################################################
    def gen_piece(self, phys, sub, piece):
    ###########################################################################
        filepath, was_filegen, insert_regex, self_begin_regex, self_end_regex \
            = [item(phys, sub, self) for item in PIECES[piece]]

        if was_filegen:
            # If freshly generated file, we're done
            pass
        else:
            orig_lines = filepath.open().readlines()
            needs_rewrite = False
            gen_lines  = getattr(self, "gen_{}".format(piece))(phys, sub).splitlines()

            # Check to see if piece already exists
            try:
                existing_piece_line_range = check_existing_piece(orig_lines, self_begin_regex, self_end_regex)
            except Exception as e:
                expect(False, "Problem parsing file {} for existing piece {}: {}".format(filepath, piece, e))

            if existing_piece_line_range is not None:
                # Replace existing
                if self._dry_run:
                    print("In file {}, would replace:\n{}\n\nWITH:\n{}".\
                          format(filepath, "\n".join(orig_lines[slice(*existing_piece_line_range)]), "\n".join(gen_lines)))
                elif not self._overwrite:
                    print("Already detected piece {} for subroutine {} in file {}, code:\n{}".\
                          format(piece, sub, filepath, "\n".join(orig_lines[slice(*existing_piece_line_range)])))
                else:
                    orig_lines[slice(*existing_piece_line_range)] = gen_lines
                    needs_rewrite = True
            else:
                # Look for place to insert and insert
                insert_line = find_insertion(orig_lines, insert_regex)
                expect(insert_line is not None,
                       "Could not find place to insert generated code for {} in file {} based on regex {}".\
                       format(piece, filepath, insert_regex))

                if self._dry_run:
                    print("In file {}, at line {}, would insert:\n{}}".format(filepath, insert_line, "\n".join(gen_lines)))
                else:
                    orig_lines[insert_line:insert_line] = gen_lines
                    needs_rewrite = True

            if needs_rewrite:
                with filepath.open("w") as fd:
                    fd.writelines(orig_lines)

    ###########################################################################
    def gen_f90_c2f_bind(self, phys, sub):
    ###########################################################################
        arg_data = self._get_arg_data(phys, sub)
        arg_names = ", ".join([item[ARG_NAME] for item in arg_data])
        arg_decls = gen_arg_f90_decls(arg_data)
        phys_mod = "micro_p3" if phys == "p3" else phys
        result = \
"""  subroutine {sub}_c({arg_names}) bind(C)
    use {phys_mod}, only : {sub}

    {arg_decls}

    call {sub}({arg_names})
  end subroutine {sub}_c
""".format(sub=sub, arg_names=arg_names, phys_mod=phys_mod, arg_decls="\n    ".join(arg_decls))

        return result

    ###########################################################################
    def gen_f90_f2c_bind(self, phys, sub):
    ###########################################################################
        arg_data = self._get_arg_data(phys, sub)
        arg_names = ", ".join([item[ARG_NAME] for item in arg_data])
        arg_decls = gen_arg_f90_decls(arg_data)
        result = \
"""  subroutine {sub}_f({arg_names}) bind(C)
    use iso_c_binding

    {arg_decls}
  end subroutine {sub}_f
""".format(sub=sub, arg_names=arg_names, arg_decls="\n    ".join(arg_decls))

        return result

    ###########################################################################
    def gen_cxx_c2f_bind_decl(self, phys, sub):
    ###########################################################################
        arg_data  = self._get_arg_data(phys, sub)
        arg_decls = gen_arg_cxx_decls(arg_data)
        result = "void {sub}_c({arg_sig});\n".format(sub=sub, arg_sig=", ".join(arg_decls))
        return result

    ###########################################################################
    def gen_cxx_c2f_glue_decl(self, phys, sub):
    ###########################################################################
        struct_name      = get_data_struct_name(sub)
        result = "void {sub}({struct_name}& d);".format(sub=sub, struct_name=struct_name)
        return result

    ###########################################################################
    def gen_cxx_c2f_glue_impl(self, phys, sub):
    ###########################################################################
        arg_data         = self._get_arg_data(phys, sub)
        arg_data_args    = gen_cxx_data_args(arg_data)
        need_transpose   = needs_transpose(arg_data)
        transpose_code_1 = "\n  d.transpose<ekat::util::TransposeDirection::c2f>();\n" if need_transpose else ""
        transpose_code_2 = "\n  d.transpose<ekat::util::TransposeDirection::f2c>();\n" if need_transpose else ""
        data_struct      = get_data_struct_name(sub)
        init_code        = get_physics_data(phys, INIT_CODE)

        result = \
"""void {sub}({data_struct}& d)
{{
  {init_code}{transpose_code_1}
  {sub}_c({arg_data_args});{transpose_code_2}
}}
""".format(sub=sub, data_struct=data_struct, init_code=init_code, transpose_code_1=transpose_code_1, transpose_code_2=transpose_code_2, arg_data_args=arg_data_args)
        return result

    ###########################################################################
    def gen_cxx_c2f_data(self, phys, sub):
    ###########################################################################
        arg_data         = self._get_arg_data(phys, sub)
        struct_members   = gen_struct_members(arg_data)
        any_arrays       = has_arrays(arg_data)
        struct_name      = get_data_struct_name(sub)
        inheritance      = " : public PhysicsTestData" if any_arrays else ""
        api              = gen_struct_api(phys, struct_name, arg_data) if any_arrays else ""

        result = \
"""struct {struct_name}{inheritance} {{
{struct_members}{api}
}};

""".format(struct_name=struct_name, inheritance=inheritance, struct_members=struct_members, api="\n  ".join(api))
        return result

    ###########################################################################
    def gen_cxx_f2c_bind_decl(self, phys, sub):
    ###########################################################################
        arg_data  = self._get_arg_data(phys, sub)
        arg_decls = gen_arg_cxx_decls(arg_data)

        return "void {sub}_f({arg_sig});".format(sub=sub, arg_sig=", ".join(arg_decls))

    ###########################################################################
    def gen_cxx_f2c_bind_impl(self, phys, sub):
    ###########################################################################
        decl = self.gen_cxx_f2c_bind_decl(phys, sub).rstrip(";")

        # TODO - is it possible to fill-out some or all of the implementation?
        result = \
"""{decl}
{{
  // TODO
}} // end {sub}_f
""".format(sub=sub, decl=decl)
        return result

    ###########################################################################
    def gen_cxx_func_decl(self, phys, sub):
    ###########################################################################
        arg_data = self._get_arg_data(phys, sub)
        arg_decls = gen_arg_cxx_decls(arg_data, kokkos=True)

        return "void {sub}({arg_sig});".format(sub=sub, arg_sig=", ".join(arg_decls))

    ###########################################################################
    def gen_cxx_func_impl(self, phys, sub):
    ###########################################################################
        decl = self.gen_cxx_func_decl(phys, sub).rstrip(";")

        # I don't think any intelligent guess at an impl is possible here
        result = \
"""{decl}
{{
  // TODO
}}
""".format(decl=decl)
        return result

    ###########################################################################
    def gen_cxx_bfb_unit_decl(self, phys, sub):
    ###########################################################################
        test_struct = get_data_test_struct_name(sub)
        return "struct {};".format(test_struct)

    ###########################################################################
    def gen_cxx_bfb_unit_impl(self, phys, sub):
    ###########################################################################
        arg_data = self._get_arg_data(phys, sub)
        data_struct = get_data_struct_name(sub)
        has_array = has_arrays(arg_data)
        need_transpose = needs_transpose(arg_data)

        gen_random = "" if not has_array else \
"""

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize();
    }"""

        c2f_transpose_code = "" if not need_transpose else \
"""
      d.transpose<ekat::util::TransposeDirection::c2f>(); // _f expects data in fortran layout"""
        f2c_transpose_code = "" if not need_transpose else \
"""
      d.transpose<ekat::util::TransposeDirection::f2c>(); // go back to C layout"""

        _, _, _, scalars, ik_reals, ij_reals, i_reals, i_ints = group_data(arg_data, filter_out_intent="in")
        check_scalars, check_arrays = "", ""
        i_arrays = i_reals + i_ints
        for scalar in scalars:
            check_scalars += "    REQUIRE(f90_data.{name} == cxx_data.{name});\n".format(name=scalar)

        for items, total_call in [(ik_reals, "total1x2()"), (ij_reals, "total1x3()"), (i_arrays, "dim1")]:
            if len(items) > 0:
                check_arrays += "      for (Int k = 0; k < f90_data.{}; ++k) {{\n".format(total_call)
                for item in items:
                    check_arrays += "        REQUIRE(f90_data.{name}[k] == cxx_data.{name}[k]);\n".format(name=item)
                check_arrays += "      }\n"

        result = \
"""  static void bfb()
  {{
    {data_struct} f90_data[] = {{
      // TODO
    }};

    static constexpr Int num_runs = sizeof(f90_data) / sizeof({data_struct});{gen_random}

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    {data_struct} cxx_data[] = {{
      // TODO
    }};

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {{
      // expects data in C layout
      calc_shoc_vertflux(d);
    }}

    // Get data from cxx
    for (auto& d : cxx_data) {{{c2f_transpose_code}
      calc_shoc_vertflux_f(d.shcol(), d.nlev(), d.nlevi(), d.tkh_zi, d.dz_zi, d.invar, d.vertflux);{f2c_transpose_code}
    }}

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {{
      {data_struct}& d_f90 = f90_data[i];
      {data_struct}& d_cxx = cxx_data[i];
{check_scalars}{check_arrays}
    }}
}}
""".format(data_struct=data_struct,
           gen_random=gen_random,
           c2f_transpose_code=c2f_transpose_code,
           f2c_transpose_code=f2c_transpose_code,
           check_scalars=check_scalars,
           check_arrays=check_arrays)

        return result

    ###########################################################################
    def gen_cxx_eti(self, phys, sub):
    ###########################################################################
        include_file = PIECES["cxx_func_impl"][0](phys, sub, self)

        result = \
"""#include {include_file}

namespace scream {{
namespace {phys} {{

/*
 * Explicit instantiation for doing {sub} on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace {phys}
} // namespace scream
""".format(include_file=include_file, phys=phys)

        return result

    ###########################################################################
    def gen_boiler(self):
    ###########################################################################
        all_success = True
        for sub in self._subs:
            for phys in self._physics:
                for piece in self._pieces:
                    try:
                        gen_piece(phys, sub, piece)
                    except Exception as e:
                        print("Warning: failed to generate subroutine {} piece {} for physics {}, error: {}".\
                              format(sub, piece, phys, e))
                        all_success = False

        return all_success
