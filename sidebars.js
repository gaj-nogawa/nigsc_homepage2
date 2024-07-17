/**
 * Creating a sidebar enables you to:
 - create an ordered group of docs
 - render a sidebar for each doc of that group
 - provide next/previous navigation
 The sidebars can be generated from the filesystem, or explicitly defined here.
 Create as many sidebars as you want.
*/

module.exports = {
    // By default, Docusaurus generates a sidebar from the docs folder structure
    //tutorialSidebar: [{type: 'autogenerated', dirName: '.'}],

    // But you can create a sidebar manually
    /*
      tutorialSidebar: [
      {
      type: 'category',
      label: 'Tutorial',
      items: ['hello'],
      },
      ],
    */

    topSidebar: [
        {
            type: "link",
            label: "お知らせ",
            href: "/blog/tags/news"
        },
        {
            type: "link",
            label: "メンテナンス情報",
            href: "/blog/tags/maintenance"
        },
        {
            type: "category",
            label: "システム構成",
            items: [
                "guides/top",
                {
                    type: "link",
                    label: "サイト内検索",
                    href: "https://sc.ddbj.nig.ac.jp/search"
                },

                {
                    type: "category",
                    label: "情報セキュリティ方針(ISO27001)",
                    link: {
                        type: "doc",
                        id: "guides/security-policy",
                    },
                    items: [
                        "guides/ISMS_Certificate",
                    ],
                },
                "guides/overview",
                "guides/hardware",
                {
                    type: "category",
                    label: "ソフトウェア",
                    link: {
                        type: "doc",
                        id: "software/software",
                    },
                    items: [
                        "software/software_update_log",
                    ],
                },
            ],
        },
        "guides/start_of_use",
        "guides/introduction",
        /*{
            type: "doc",
            id: "guides/start_of_use",
            label: "アカウントの新規登録から利用開始までの流れ",
        },
        {
            type: "doc",
            id: "guides/introduction",
            label: "ログイン方法の概要",
        },*/
        {
            type: "category",
            label: "一般解析区画の使い方",
            items: [
                "general_analysis_division/ga_introduction",
                "general_analysis_division/ga_application",
                "general_analysis_division/ga_login",
                "general_analysis_division/ga_usage",
                "general_analysis_division/ga_queue",
                "general_analysis_division/ga_transfer",
                "general_analysis_division/advance_reservation",
                "general_analysis_division/ga_lustre",
                "general_analysis_division/largescale_storage",
                "operation/qstatGC",
                "operation/gfree",
            ],
        },
        {
            type: "category",
            label: "個人ゲノム解析区画の使い方",
            items: [
                "personal_genome_division/pg_application",
                {
                    type: "category",
                    label: "ログイン方法 (個人ゲノム解析区画)",
                    link: {
                        type: "doc",
                        id: "personal_genome_division/pg_login",
                    },
                    items: [
                        {
                            type: "category",
                            label: "SSL-VPNクライアントソフトウェアのインストール方法",
                            items: [
                                "personal_genome_division/pg_login_ssl-vpn_install_win",
                                "personal_genome_division/pg_login_ssl-vpn_install_mac",
                                "personal_genome_division/pg_login_ssl-vpn_install_linux",
                            ],
                        },
                        {
                            type: "category",
                            label: "SSL-VPNクライアントの設定方法",
                            items: [
                                "personal_genome_division/pg_login_ssl-vpn_configure_file_win",
                                "personal_genome_division/pg_login_ssl-vpn_configure_file_mac",
                                "personal_genome_division/pg_login_ssl-vpn_configure_file_linux",
                            ],
                        },
                        {
                            type: "category",
                            label: "SSL-VPNへの接続方法",
                            items: [
                                "personal_genome_division/pg_login_ssl-vpn_connection_win",
                                "personal_genome_division/pg_login_ssl-vpn_connection_mac",
                                "personal_genome_division/pg_login_ssl-vpn_connection_linux",
                            ],
                        },
                    ],
                },
                "personal_genome_division/pg_usage",
                "personal_genome_division/gpu_slurm",
                "personal_genome_division/pg_transfer",
                "general_analysis_division/ga_lustre",
                "general_analysis_division/largescale_storage",
                "personal_genome_division/group_cloud",
            ],
        },
    ],


    softwareSidebar: [
        {
            type: 'category',
            label: "ジョブスケジューラ",
            items: [
                {
                    type: 'category',
                    label: 'Grid Engine',
                    items: [
                        "software/grid_engine/grid_engine",
                        "software/grid_engine/interactive_jobs",
                        "software/grid_engine/batch_jobs",
                        "software/grid_engine/parallel_jobs/parallel_jobs",
                        "software/grid_engine/array_jobs",
                        "software/grid_engine/other_commands",
                        "software/qsub_beta",
                    ]
                },
                {
                    type: 'category',
                    label: 'Slurm',
                    items: [
                        "software/Slurm/Slurm",
                        "software/Slurm/batch_jobs",
                        "software/Slurm/interactive_jobs",
                        "software/Slurm/parallel_jobs",
                        "software/Slurm/array_jobs",
                        "software/Slurm/other_commands",
                    ]
                },


            ]
        },
        {
            type: 'category',
            label: 'パッケージマネージャ',
            items: [
                {
                    type: "category",
                    label: "spack",
                    items: [
                        "software/spack/install_spack",
                        "software/spack/use_spack",
                    ],
                },
                {
                    type: "category",
                    label: "Environmental Modules",
                    items: [
                        "software/environmental_modules/environmental_modules",
                    ],
                },

            ]
        },


        {
            type: 'category',
            label: 'コンテナ・解析パイプライン',
            items: [
                "software/BioContainers/BioContainers",
                "software/Apptainer/Apptainer",
            ]
        },

        {
            type: 'category',
            label: "データ転送・データ共有",
            items: [
                {
                    type: 'category',
                    label: 'Aspera',
                    items: [
                        "software/aspera/aspera",
                        "software/aspera/install_Aspera",
                    ],
                },
                {
                    type: "category",
                    label: "HCP tools",
                    items: [
                        "software/Archaea_tools/Archaea_tools",
                        "software/Archaea_tools/hcptools_conf",
                    ],
                }
            ]
        },

        {
            type: 'category',
            label: "開発環境・ライブラリ",
            items: [
                "software/python",
                "software/R/R",
                "software/R/r_studio_server",
                "software/jupyter_notebook",
                "software/jupyter_lab",
                "software/java",
                "software/typescript",
                "software/rust",
                "software/gcc/gcc",
                {
                    type: "category",
                    label: "C/C++の使い方 (Intel Compiler)",
                    link: {
                        type: "doc",
                        id: "software/intel_compiler/intel_compiler",
                    },
                    items: [
                        "software/intel_compiler/intel_compiler_ParallelComLib",
                        "software/intel_compiler/intel_compiler_NumCalc",
                        "software/intel_compiler/intel_compiler_AILib",
                        "software/intel_compiler/intel_compiler_profiler",
                        "software/intel_compiler/intel_compiler_debugger",
                        "software/intel_compiler/intel_compiler_OtherLanguages",
                    ],
                },
                "software/pgi_compiler/pgi_compiler",
                "software/cuda",
                "software/go",
            ]

        }
    ],


    applicationsSidebar : [
        {
            type: "category",
            label: "利用規定等",
            items: [
                "application/application",
                {
                    type: "category",
                    label: "誓約書に署名する方法",
                    items: [
                        "application/signing_PDF",
                        "application/signing_PDF_domestic_resident",
                        "application/signing_PDF_non-resident",
                    ]
                },
                "application/use_policy",
                "application/legislation",
            ],
        },
        "application/registration",
        "application/change_account_info",
        "application/renewal",
        "application/resource_extension",
        {
            type: "category",
            label: "パスワード・公開鍵の設定方法",
            items: [
                "application/ssh_keys",
                {
                    type: "category",
                    label: "SSH公開鍵・秘密鍵の生成方法",
                    items: [
                        "application/ssh_keys_ssh-keygen_win",
                        "application/ssh_keys_ssh-keygen_mac",
                        "application/ssh_keys_ssh-keygen_linux",
                    ],
                },
                {
                    type: "category",
                    label: "遺伝研スパコンへのSSH公開鍵の登録・変更方法",
                    items: [
                        "application/ssh_keys_register_win",
                        "application/ssh_keys_register_mac",
                        "application/ssh_keys_register_linux",
                    ],
                },
                "application/ssh_keys_connect_NIGsupercomputer",
                "application/change_loginpwd",
           ],
        },
        {
           type: "category",
           label: "課金サービス利用方法",
           items: [
            "application/billing_service",
            "application/resource_extension",
            "application/invoice",
            ],
        },
        {
            type: "link",
            label: "よくある質問(FAQ)",
            href: "/faq/faq_software"
         },
        {
            type: "link",
            label: "Github Discussions(Q&A)",
            href: "https://github.com/nig-sc/nigsc_homepage2/discussions"
        },
        "application/reference",
        {
        type: "link",
        label: "古いドキュメント",
        href: "/oldDocuments/software/R",
        },
    ],


    advancedGuidesSidebar : [
        {
            type: "category",
            label: "活用方法",
            items: [
                "advanced_guides/advanced_guide_2023",
                "advanced_guides/advanced_guide_2020-2022",
            ],
        },
        {
            type: "category",
            label: "Rhelixa RNAseq",
            items: [
                "advanced_guides/Rhelixa_RNAseq",
                "advanced_guides/Rhelixa_RNAseq_manual",
            ],
        },
        {
            type: "category",
            label: "Rhelixa Graphing Tool",
            items: [
                "advanced_guides/Rhelixa_RNAseq_Visualization",
            ],
        },
        {
            type: "category",
            label: "Alphafold",
            items: [
                "advanced_guides/Alphafold_2_1",
                "advanced_guides/Alphafold_2_2",
                "advanced_guides/Alphafold_2_3",
            ],
        },
        {
            type: "category",
            label: "TogoImputation (beta)",
            items: [
                "advanced_guides/imputation_server",
                "advanced_guides/imputation_server_install",
                "advanced_guides/imputation_server_tutorial",
                "advanced_guides/imputation_server_tutorial2",
                "advanced_guides/imputation_server_hail",
                "advanced_guides/imputation_server_hibag",
            ],
        },
        {
            type: "category",
            label: "NVIDIA Parabricks",
            items: [
                "advanced_guides/parabricks/parabricks",
            ],
        },
        {
            type: "category",
            label: "Benchmark",
            items: [
                "advanced_guides/benchmark_dorado",
                "advanced_guides/benchmark_parabricks",
            ],  
        },
        {
            type: "doc",
            id: "advanced_guides/experimental",
            label: "Experimental",
        },

        {
            type: "category",
            label: "講習会",
            items: [
                "advanced_guides/IIBMP2021",
                "advanced_guides/IIBMP2020",
            ],
        },
        {
            type: "category",
            label: "利用方法解説",
            items: [
                "advanced_guides/commentary",
            ],
        },

    ],

    operationInfoSidebar: [
        {
            type: "category",
            label: "稼働状況",
            items: [
                "operation/operation",
                "operation/Total_PowerConsumption",
            ]
        },
    ],

    reportSidebar: [
        "report/report",
        {
            type: "category",
            label: "論文リスト",
            items: [
            "report/papers_2022",
            "report/papers_2021",
            "report/papers_2020",
            "report/papers_2019",
            "report/papers_2018",
            "report/papers_2017",
            "report/papers_2016",
            "report/papers_2015",
            "report/papers_2014",
            "report/papers_2013",
            "report/papers_2012",
            ],
        },
        {
            type: "category",
            label: "アカウントの所属機関数",
            items: [
            "report/number_of_account_institutions",
            ]
        },
    ],


    faqSidebar: [
        {
            type: 'category',
            label: "FAQ : software",
            items: [
                "faq/faq_software",
                "faq/faq_hcptools",
                "faq/faq_aspera",
                "faq/faq_grid_engine",
            ]
        },
        {
            type: 'category',
            label: "FAQ : 一般解析区画",
            items: [
                "faq/faq_login_general",
            ]
        },
        {
            type: 'category',
            label: "FAQ : 個人ゲノム解析区画",
            items: [
                "faq/faq_login_personal",
            ]
        },
        {
            type: 'category',
            label: "FAQ : 利用規定等",
            items: [
                {
                type: 'category',
                label: "FAQ : 誓約書に署名する方法",
                items: [
                    "faq/faq_signing_PDF",
                ]
                },
            ]
        },
        {
            type: 'category',
            label: "FAQ : 各種申請",
            items: [
                "faq/faq_NewUser_registration",
                "faq/faq_renewal",
            ]
        },
        {
            type: 'category',
            label: "FAQ : パスワード・公開鍵の設定方法",
            items: [
                "faq/faq_sshkeys_windows",
            ]
        },
        {
            type: 'category',
            label: "FAQ : 課金サービス利用方法",
            items: [
                "faq/faq_billing_service",
                {
                    type: 'category',
                    label: "利用計画表の提出",
                    items: [
                        "faq/faq_change_StorageCapacity",
                    ]
                },
            ]
        },
        {
            type: "category",
            label: "FAQ：OS移行に伴うご質問",
            link: {
                type: "doc",
                id: "faq/faq_os_migration",
            },
            items: [
                "faq/faq_os_migration_login",
                "faq/faq_os_migration_env-var",
            ],
        },
        {
            type: "link",
            label: "Github Discussions(Q&A)",
            href: "https://github.com/nig-sc/nigsc_homepage2/discussions"
        },

    ],
    oldDocumentsSidebar: [
        {
            type: "category",
            label: "古いドキュメント",
            items: [
                {
                    type: 'category',
                    label: "software",
                    items: [
                        "oldDocuments/software/guix/guix",
                        "oldDocuments/software/R",
                    ],
                },
            ],
        },
        "faq/old_document",
    ],


};
