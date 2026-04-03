add_rules("mode.debug", "mode.release", "mode.releasedbg")
set_defaultmode("releasedbg")

add_requires("emhash", "thread-pool")

function setup_target(target_name)
	set_kind("binary")
	set_languages("cxx20")
	set_rundir("$(projectdir)/")
	
	-- 平台特定的编译标志
	if is_plat("windows") then
		set_toolchains("llvm")
		add_cxflags("-march=native")
	elseif is_plat("macosx", "linux") then
		set_toolchains("clang")
		add_cxflags("-march=native")
	end
	add_defines(
		"__ORDER_LITTLE_ENDIAN__=1234",
		"__ORDER_BIG_ENDIAN__=4321",
		"__BYTE_ORDER__=1234"
	)
end

target("test_SA")
	setup_target("test_SA")
	add_files("src/optimizer/sa/test_sa.cpp")
	add_includedirs("src/include")
	add_packages("emhash")

target("test_TABU")
	setup_target("test_TABU")
	add_files("src/optimizer/tabu/test_tabu.cpp")
	add_includedirs("src/include")
	add_packages("emhash")

target("test_SA_MT")
	setup_target("test_SA_MT")
	add_files("src/optimizer/sa_mt/test_sa_mt.cpp")
	add_includedirs("src/include")
	add_packages("emhash", "thread-pool")

target("test_ALNS")
	setup_target("test_ALNS")
	add_files("src/optimizer/alns/test_alns.cpp")
	add_includedirs("src/include")
	add_packages("emhash", "thread-pool")

target("test_AMHB")
	setup_target("test_AMHB")
	add_files("src/optimizer/amhb/test_amhb.cpp")
	add_includedirs("src/include")
	add_packages("emhash", "thread-pool")

