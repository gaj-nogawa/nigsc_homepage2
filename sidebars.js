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
                "guides/overview",
                "guides/security-policy",
                "guides/hardware",
                "software/software"
            ]
        },
        "start_the_service",
        "introduction",
        {
            type: "category",
            label: "一般解析区画の使い方",
            items: [
                "general_analysis_division/ga_introduction",
                "general_analysis_division/ga_application",
                "general_analysis_division/ga_login",
                "general_analysis_division/ga_usage",
                "general_analysis_division/ga_transfer",
                "general_analysis_division/advance_reservation",
                "general_analysis_division/ga_lustre",
                "general_analysis_division/largescale_storage",
            ]
        },
        {
            type: "category",
            label: "個人ゲノム解析区画の使い方",
            items: [
                "personal_genome_division/pg_introduction",
                "personal_genome_division/pg_application",
                "personal_genome_division/pg_login",
                "personal_genome_division/pg_usage",
                "personal_genome_division/pg_transfer",
                "general_analysis_division/ga_lustre",
                "general_analysis_division/largescale_storage",
                "personal_genome_division/group_cloud",
            ]
        }
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
                        "software/slurm"
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
                    label: "Guix",
                    items: [
                        "software/guix/guix",
                    ],
                },
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
                "software/R",
                "software/jupyter_notebook",
                "software/jupyter_lab",
                "software/java",
                "software/typescript",
                "software/rust",
                "software/cuda",
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
        {
            type: "category",
            label: "利用申請等",
            items: [
            "application/registration",
            "application/resource_extension",
            "application/renewal",
            ],
        },
        {
            type: "category",
            label: "パスワード・公開鍵の設定方法",
            items: [
            "application/ssh_keys",
            "application/ssh_keys_mac",
            "application/ssh_keys_windows",
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
        "faq/old_document",
    ],


    advancedGuidesSidebar : {
        "活用方法": [
            "advanced_guides/advanced_guide_2023",
            "advanced_guides/advanced_guide_2020-2022",
        ],
        "Rhelixa RNAseq": [
            "advanced_guides/Rhelixa_RNAseq",
            "advanced_guides/Rhelixa_RNAseq_manual",
        ],
        "Alphafold":[
            "advanced_guides/Alphafold_2_1",
            "advanced_guides/Alphafold_2_2",
            "advanced_guides/Alphafold_2_3",
        ],
        "NBDC-DDBJ Imputation Server (beta)": [
            "advanced_guides/imputation_server",
            "advanced_guides/imputation_server_install",
            "advanced_guides/imputation_server_tutorial",
            "advanced_guides/imputation_server_tutorial2",
            "advanced_guides/imputation_server_hail",
        ],
        "NVIDIA Parabricks": [
            "advanced_guides/parabricks/parabricks",
        ],
        "講習会": [
            "advanced_guides/IIBMP2021",
            "advanced_guides/IIBMP2020",
        ],
        "利用方法解説": [
            "advanced_guides/commentary"
        ],
    },
    operationInfoSidebar: [
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
            label: "稼働状況",
            items: [
                "operation/operation",
                "operation/qstatGC",
                "operation/gfree",
                "operation/Total_PowerConsumption",
            ]
        }
    ],
    reportSidebar: {
        "各種統計": [
            "report/report",
        ],
        "論文リスト": [
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
            label: "FAQ : 利用申請等",
            items: [
                "faq/faq_NewUser_registration",
                "faq/faq_renewal",
            ]
        },
        {
            type: 'category',
            label: "FAQ : パスワード・公開鍵の設定方法",
            items: [
                "faq/faq_sshkeys_mac",
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
    ],
};
